"""Bounded, non-interactive network enrichment for the command-line runner.

Provider HTTP work lives in a spawned process so deadlines and cancellation also
stop provider polling, retry sleeps, and all of its request threads.
"""

from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
import json
import logging
import multiprocessing
import os
import time

from .cli_config import CliError, NETWORK_API_KEY_ENV, NETWORK_SOURCES


SOURCES = NETWORK_SOURCES
PROPERTIES = ("smiles", "cas", "iupac", "mw", "name")
IDENTIFIERS = ("smiles", "cas", "name")
API_KEY_ENV = NETWORK_API_KEY_ENV
_SUPPORTED_PROPERTIES = {
    "pubchem": set(PROPERTIES), "comptox": set(PROPERTIES),
    "cactus": set(PROPERTIES), "chemspider": {"smiles", "mw", "name"},
    "chembl": {"smiles", "iupac", "mw", "name"},
    "cas": {"smiles", "cas", "mw", "name"},
}
_SUPPORTED_IDENTIFIERS = {
    source: set(IDENTIFIERS) if source != "chembl" else {"smiles", "name"}
    for source in SOURCES
}
_DEFAULTS = {
    "sources": ["pubchem"], "send": ["smiles"], "max_workers": 4,
    "timeout": 15, "deadline": 120, "retries": 2,
    "retry_delay": 1, "interval": 0.3,
}
_INFRASTRUCTURE_ERRORS = {
    "AUTHENTICATION_FAILED", "RATE_LIMITED", "REMOTE_SERVER_ERROR", "TIMEOUT",
    "CONNECTION_FAILED", "NETWORK_ERROR", "PROVIDER_ERROR", "INVALID_RESPONSE",
}


def preflight_enrich(params, roles):
    """Check capabilities and configured identifier roles without any HTTP work."""
    options = {**_DEFAULTS, **params}
    sources, send, requested = options["sources"], options["send"], options.get("request", [])
    for field, values, allowed in (("sources", sources, SOURCES),
                                   ("send", send, IDENTIFIERS),
                                   ("request", requested, PROPERTIES)):
        if not isinstance(values, list) or not values or any(value not in allowed for value in values):
            raise CliError("INVALID_ENRICH_PARAMETERS", f"enrich.{field} requires supported values.",
                           {"field": field, "allowed": list(allowed)})
    for identifier in send:
        if not roles.get(identifier):
            raise CliError("COLUMN_REQUIRED", f"enrich requires a column for '{identifier}'.",
                           {"role": identifier})
        if not any(identifier in _SUPPORTED_IDENTIFIERS[source] for source in sources):
            raise CliError("UNSUPPORTED_IDENTIFIER", "The selected sources cannot query an identifier type.",
                           {"identifier_type": identifier, "sources": sources})
    for prop in requested:
        if not any(prop in _SUPPORTED_PROPERTIES[source]
                   and set(send) & _SUPPORTED_IDENTIFIERS[source] for source in sources):
            raise CliError("UNSUPPORTED_PROPERTY", "The selected sources cannot supply a requested property.",
                           {"property": prop, "sources": sources, "send": send})
    missing = [API_KEY_ENV[source] for source in sources
               if source in API_KEY_ENV and not os.environ.get(API_KEY_ENV[source], "").strip()]
    if missing:
        raise CliError("MISSING_API_KEY", "Selected enrichment sources require environment credentials.",
                       {"environment_variables": missing})


def _error(code, source, retryable=False, http_status=None):
    return {"code": code, "source": source, "retryable": retryable, "http_status": http_status}


class _SafeRequestError(Exception):
    def __init__(self, error):
        super().__init__(error["code"])
        self.error = error


def _classify_exception(exc, source):
    import requests

    if isinstance(exc, _SafeRequestError):
        return exc.error
    if isinstance(exc, requests.exceptions.Timeout):
        return _error("TIMEOUT", source, True)
    if isinstance(exc, requests.exceptions.ConnectionError):
        return _error("CONNECTION_FAILED", source, True)
    if isinstance(exc, requests.exceptions.HTTPError) and exc.response is not None:
        return _classify_status(exc.response.status_code, source)
    if isinstance(exc, requests.exceptions.RequestException):
        return _error("NETWORK_ERROR", source, True)
    if isinstance(exc, (ValueError, TypeError, KeyError, IndexError)):
        return _error("INVALID_RESPONSE", source)
    return _error("PROVIDER_ERROR", source)


def _classify_status(status, source):
    if status in (401, 403):
        code, retryable = "AUTHENTICATION_FAILED", False
    elif status == 404:
        code, retryable = "NOT_FOUND", False
    elif status == 400:
        code, retryable = "INVALID_IDENTIFIER", False
    elif status == 429:
        code, retryable = "RATE_LIMITED", True
    elif status >= 500:
        code, retryable = "REMOTE_SERVER_ERROR", True
    else:
        code, retryable = "NETWORK_ERROR", False
    return _error(code, source, retryable, int(status))


class _RequestPacer:
    def __init__(self, interval):
        from threading import Lock

        self.lock = Lock()
        self.auth_lock = Lock()
        self.auth_failures = {}
        self.interval = interval
        self.next_request = 0.0

    def wait(self):
        with self.lock:
            delay = self.next_request - time.monotonic()
            if delay > 0:
                time.sleep(delay)
            self.next_request = time.monotonic() + self.interval

    def authentication_failure(self, source):
        with self.auth_lock:
            return self.auth_failures.get(source)

    def remember_authentication_failure(self, error):
        with self.auth_lock:
            self.auth_failures[error["source"]] = error


def _configure_transport(service, options, errors, source_state, pacer):
    """Override legacy fixed timeouts and adapter retries in this worker only."""
    from requests.adapters import HTTPAdapter

    for adapter in service.session.adapters.values():
        adapter.close()
    service.session.mount("https://", HTTPAdapter(max_retries=0))
    service.session.mount("http://", HTTPAdapter(max_retries=0))
    original_request = service.session.request

    def request(method, url, **kwargs):
        kwargs["timeout"] = options["timeout"]
        for attempt in range(options["retries"] + 1):
            pacer.wait()
            cached_failure = pacer.authentication_failure(source_state[0])
            if cached_failure:
                errors.append(cached_failure)
                raise _SafeRequestError(cached_failure)
            response = None
            try:
                response = original_request(method, url, **kwargs)
                if response.status_code < 400:
                    return response
                failure = _classify_status(response.status_code, source_state[0])
            except Exception as exc:
                failure = _classify_exception(exc, source_state[0])
            if failure["retryable"] and attempt < options["retries"]:
                if response is not None:
                    response.close()
                time.sleep(options["retry_delay"])
                continue
            errors.append(failure)
            if failure["code"] == "AUTHENTICATION_FAILED":
                pacer.remember_authentication_failure(failure)
            if response is not None and response.status_code == 404:
                # Some providers have a legitimate second lookup after a 404.
                return response
            if response is not None:
                response.close()
            raise _SafeRequestError(failure) from None

    service.session.request = request


def _clean_identifier(value, identifier_type, service):
    if value is None or not str(value).strip():
        return None
    value = str(value).strip()
    if identifier_type == "cas":
        return service._validate_and_clean_cas(value)
    if identifier_type == "smiles":
        from rdkit import Chem

        return value if Chem.MolFromSmiles(value) is not None else None
    return value


def _enrich_row(identifiers, options, pacer):
    from .web_request import WebService

    service = WebService(
        sources=options["sources"], force_api_mode=True, batch_limit=0,
        interval=0, retries=0, delay=0,
        **{f"{source}_api_key": os.environ.get(variable)
           for source, variable in API_KEY_ENV.items()},
    )
    values, provenance, errors = {}, {}, []
    source_state = [None]
    _configure_transport(service, options, errors, source_state, pacer)
    valid_identifier = False
    disabled_sources = set()
    try:
        for identifier_type in options["send"]:
            identifier = _clean_identifier(identifiers.get(identifier_type), identifier_type, service)
            if identifier is None:
                continue
            valid_identifier = True
            for source in options["sources"]:
                if source in disabled_sources or identifier_type not in _SUPPORTED_IDENTIFIERS[source]:
                    continue
                cached_failure = pacer.authentication_failure(source)
                if cached_failure:
                    errors.append(cached_failure)
                    disabled_sources.add(source)
                    continue
                source_state[0] = source
                # Fetch separately so a later endpoint failure cannot discard a
                # property already supplied successfully by this same provider.
                for prop in options["request"]:
                    if prop in values or prop not in _SUPPORTED_PROPERTIES[source]:
                        continue
                    try:
                        result = service._fetch_functions[source](identifier, {prop}, identifier_type) or {}
                        value = result.get(prop)
                        if value is not None and str(value).strip():
                            values[prop] = str(value)
                            provenance[prop] = {"source": source, "identifier_type": identifier_type}
                    except Exception as exc:
                        failure = _classify_exception(exc, source)
                        errors.append(failure)
                        if failure["code"] == "AUTHENTICATION_FAILED":
                            pacer.remember_authentication_failure(failure)
                            disabled_sources.add(source)
                            break
                        if failure["code"] == "INVALID_IDENTIFIER":
                            break
                if len(values) == len(options["request"]):
                    break
            if len(values) == len(options["request"]):
                break
    finally:
        service.session.close()
    unique_errors = list({(error["code"], error["source"], error["http_status"]): error
                          for error in errors}.values())
    missing = [prop for prop in options["request"] if prop not in values]
    if not missing:
        status = "complete"
    elif values:
        status = "partial"
    elif not valid_identifier or (unique_errors and all(
            error["code"] == "INVALID_IDENTIFIER" for error in unique_errors)):
        status = "invalid_identifier"
    elif any(error["code"] in _INFRASTRUCTURE_ERRORS for error in unique_errors):
        status = "failed"
    else:
        status = "not_found"
    return {"values": values, "provenance": provenance, "errors": unique_errors,
            "missing": missing, "status": status}


def _network_worker(connection, rows, options):
    # Provider and SDK diagnostics can include credentials or complete URLs.
    # Only explicitly structured results cross this process boundary.
    with open(os.devnull, "w", encoding="utf-8") as sink, redirect_stdout(sink), redirect_stderr(sink):
        logging.disable(logging.CRITICAL)
        try:
            from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
            from rdkit import RDLogger

            RDLogger.DisableLog("rdApp.*")
            pacer = _RequestPacer(options["interval"])
            pending_rows = iter(enumerate(rows))
            with ThreadPoolExecutor(max_workers=options["max_workers"]) as executor:
                active = {}

                def submit_next():
                    item = next(pending_rows, None)
                    if item is not None:
                        index, row = item
                        active[executor.submit(_enrich_row, row, options, pacer)] = index

                for _ in range(options["max_workers"]):
                    submit_next()
                while active:
                    done, _ = wait(active, return_when=FIRST_COMPLETED)
                    for future in done:
                        index = active.pop(future)
                        connection.send(("row", index, future.result()))
                        submit_next()
            connection.send(("done",))
        except BaseException:
            try:
                connection.send(("error",))
            except (BrokenPipeError, EOFError, OSError):
                pass
        finally:
            connection.close()


def _collect_results(rows, options, *, worker=_network_worker):
    context = multiprocessing.get_context("spawn")
    receiving, sending = context.Pipe(duplex=False)
    process = context.Process(target=worker, args=(sending, rows, options), daemon=True)
    started = False
    end_time = time.monotonic() + options["deadline"]
    results = [None] * len(rows)
    try:
        process.start()
        started = True
        sending.close()
        while True:
            remaining = end_time - time.monotonic()
            if remaining <= 0:
                raise CliError("DEADLINE_EXCEEDED", "The enrichment step exceeded its total deadline.",
                               {"deadline": options["deadline"]}, retryable=True, exit_code=4)
            if receiving.poll(min(0.05, remaining)):
                try:
                    message = receiving.recv()
                except (EOFError, OSError):
                    raise CliError("ENRICH_WORKER_FAILED", "The enrichment worker stopped unexpectedly.",
                                   exit_code=4) from None
                if message[0] == "row":
                    results[message[1]] = message[2]
                elif message[0] == "done" and all(result is not None for result in results):
                    return results
                else:
                    raise CliError("ENRICH_WORKER_FAILED", "The enrichment worker could not complete the request.",
                                   exit_code=4)
            elif not process.is_alive():
                raise CliError("ENRICH_WORKER_FAILED", "The enrichment worker stopped unexpectedly.",
                               exit_code=4)
    finally:
        # This also runs on KeyboardInterrupt; no request thread survives in
        # the CLI process and no worker can continue after cancellation.
        if started:
            if process.is_alive():
                process.terminate()
            process.join(timeout=1)
            if process.is_alive():
                process.kill()
                process.join(timeout=1)
            if not process.is_alive():
                process.close()
        receiving.close()
        sending.close()


def enrich_pipeline(pipeline, params):
    """Enrich a loaded pipeline and return JSON-safe aggregate field/row counts."""
    options = {**_DEFAULTS, **params}
    options["request"] = list(dict.fromkeys(options["request"]))
    options["sources"] = list(dict.fromkeys(options["sources"]))
    options["send"] = list(dict.fromkeys(options["send"]))
    roles = {role: getattr(pipeline, f"{role}_col", None) for role in IDENTIFIERS}
    if pipeline._preprocess_key:
        roles["smiles"] = "Canonical SMILES"
    preflight_enrich(options, roles)
    for identifier in options["send"]:
        if roles[identifier] not in pipeline.df.columns:
            raise CliError("COLUMN_NOT_FOUND", "An enrichment identifier column does not exist.",
                           {"role": identifier, "column": roles[identifier]})

    import pandas as pd

    rows = [{identifier: None if pd.isna(row[roles[identifier]]) else str(row[roles[identifier]])
             for identifier in options["send"]} for _, row in pipeline.df.iterrows()]
    results = _collect_results(rows, options) if rows else []
    statuses = Counter(result["status"] for result in results)
    error_counts = Counter(error["code"] for result in results for error in result["errors"])
    requested = len(rows) * len(options["request"])
    resolved = sum(len(result["values"]) for result in results)
    stats = {
        "requested_fields": requested, "resolved_fields": resolved,
        "missing_fields": requested - resolved,
        "complete_rows": statuses["complete"], "partial_rows": statuses["partial"],
        "failed_rows": len(rows) - statuses["complete"] - statuses["partial"],
        "error_counts": dict(error_counts),
    }
    if rows and not resolved and _INFRASTRUCTURE_ERRORS.intersection(error_counts):
        raise CliError("ENRICH_FAILED", "No requested fields were resolved and network/provider failures occurred.",
                       stats, retryable=any(error["retryable"] for result in results
                                            for error in result["errors"]), exit_code=4)

    pipeline._save_checkpoint()
    before = pipeline.df.copy()
    for prop in options["request"]:
        column = f"{prop}_from_web"
        pipeline.df[column] = pd.array([result["values"].get(prop) for result in results], dtype="string")
    metadata = {
        "Query_Status": [result["status"] for result in results],
        "Query_Missing_Fields": [json.dumps(result["missing"]) for result in results],
        "Query_Errors": [json.dumps(result["errors"]) for result in results],
        "Data_Source": [json.dumps(result["provenance"]) for result in results],
        "Query_Method": [json.dumps(list(dict.fromkeys(
            value["identifier_type"] for value in result["provenance"].values()))) for result in results],
    }
    for column, values in metadata.items():
        pipeline.df[column] = pd.array(values, dtype="string")
    for role in IDENTIFIERS:
        if role in options["request"] and not (role == "smiles" and pipeline._preprocess_key):
            setattr(pipeline, f"{role}_col", f"{role}_from_web")
    pipeline._record_step("Enrich", before, pipeline.df,
                          f"Resolved fields: {resolved}/{requested}; complete rows: {statuses['complete']}/{len(rows)}")
    return stats
