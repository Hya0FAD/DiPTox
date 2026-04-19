# diptox/web_request.py
import requests
from typing import Optional, Tuple, List, Dict, Any, Set, Callable, Union
import time
import re
from threading import Lock
from tqdm import tqdm
from rdkit import Chem
from rdkit import RDLogger
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import quote
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import logging
from .logger import log_manager
logger = log_manager.get_logger(__name__)


try:
    from packaging.version import parse as parse_version
except ImportError:
    def parse_version(version_str):
        if not version_str:
            return (0, 0, 0)
        clean_str = re.sub(r'[^0-9.]', '', version_str)
        return tuple(map(int, clean_str.split('.'))) if clean_str else (0, 0, 0)

try:
    import pubchempy
    _PUBCHEMPY_AVAILABLE = True
    logging.getLogger('pubchempy').setLevel(logging.WARNING)
except ImportError:
    _PUBCHEMPY_AVAILABLE = False

try:
    from chemspipy import ChemSpider
    _CHEMSPYPY_AVAILABLE = True
except ImportError:
    _CHEMSPYPY_AVAILABLE = False

try:
    import ctxpy as ctx
    from importlib.metadata import version

    _current_ctx_version = version("ctx-python")

    if parse_version(_current_ctx_version) >= parse_version("0.0.1a10"):
        _CTXPY_AVAILABLE = True
    else:
        _CTXPY_AVAILABLE = False
        logger.warning(
            f"Detected 'ctx-python' version {_current_ctx_version}, "
            f"which is older than the required 0.0.1a10. CompTox SDK mode disabled."
        )
except ImportError:
    _CTXPY_AVAILABLE = False


class DataSource:
    PUBCHEM = 'pubchem'
    CHEMSPIDER = 'chemspider'
    COMPTOX = 'comptox'
    CACTUS = 'cactus'
    CHEMBL = 'chembl'
    CAS = 'cas'


class WebService:
    """Handles all network request operations."""

    def __init__(self, sources: Union[str, List[str]] = DataSource.PUBCHEM,
                 interval: int = 0.3,
                 retries: int = 3,
                 delay: int = 30,
                 max_workers: int = 4,
                 batch_limit: int = 1500,
                 rest_duration: int = 300,
                 chemspider_api_key: Optional[str] = None,
                 comptox_api_key: Optional[str] = None,
                 cas_api_key: Optional[str] = None,
                 force_api_mode: bool = False,
                 status_callback: Optional[Callable[[str], None]] = None):
        """
        Initializes the WebService class
        :param sources: Data source interface
        :param interval: Time interval in seconds
        :param retries: Number of retry attempts on failure.
        :param delay: Delay between retries (in seconds).
        :param max_workers: Maximum number of concurrent requests.
        :param batch_limit: Number of requests before taking a break.
        :param rest_duration: Duration of the break in seconds.
        :param chemspider_api_key: Chemspider API key.
        :param comptox_api_key: Comptox API key.
        :param cas_api_key: CAS API key.
        :param force_api_mode: Forces API mode.
        :param status_callback: Callback function for status checking.
        """
        if isinstance(sources, str):
            self.sources = [sources.lower()]
        else:
            self.sources = [source.lower() for source in sources]
        self.interval = interval
        self.retries = retries
        self.delay = delay
        self.max_workers = max_workers
        self.batch_limit = batch_limit
        self.rest_duration = rest_duration
        self.request_count = 0
        self._counter_lock = Lock()
        self.status_callback = status_callback
        self._chemspider_api_key = chemspider_api_key
        self._comptox_api_key = comptox_api_key
        self._cas_api_key = cas_api_key

        self.session = requests.Session()
        retry_strategy = Retry(
            total=3,
            connect=3,
            read=3,
            status=0,
            backoff_factor=0.1,
            allowed_methods=["HEAD", "GET", "OPTIONS", "POST"]
        )

        adapter = HTTPAdapter(
            max_retries=retry_strategy,
            pool_connections=100,
            pool_maxsize=100
        )

        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        })

        self._fetch_functions = {}
        for source in self.sources:
            source_selector = {
                DataSource.PUBCHEM: (_PUBCHEMPY_AVAILABLE, self._fetch_via_pubchem_sdk, self._fetch_via_pubchem_api),
                DataSource.CHEMSPIDER: (
                    _CHEMSPYPY_AVAILABLE, self._fetch_via_chemspider_sdk, self._fetch_via_chemspider_api),
                DataSource.COMPTOX: (_CTXPY_AVAILABLE, self._fetch_via_comptox_sdk, self._fetch_via_comptox_api),
                DataSource.CACTUS: (False, None, self._fetch_via_cactus_api),
                DataSource.CHEMBL: (False, None, self._fetch_via_chembl_api),
                DataSource.CAS: (False, None, self._fetch_via_cas_api)
            }
            if source in source_selector:
                is_sdk_available, sdk_func, api_func = source_selector[source]
                if is_sdk_available and sdk_func and not force_api_mode:
                    self._fetch_functions[source] = sdk_func
                    logger.info(f"'{source}' SDK has been found. The SDK mode will be used for operation.")
                else:
                    self._fetch_functions[source] = api_func
                    logger.info(f"'{source}' SDK not found. The program will run using the direct API mode.")
            else:
                logger.warning(f"An invalid data source has been specified: {source}")
                continue

        self._supported_props = {
            DataSource.PUBCHEM: {'smiles', 'cas', 'iupac', 'mw', 'name'},
            DataSource.CHEMSPIDER: {'smiles', 'mw', 'name'},
            DataSource.COMPTOX: {'smiles', 'cas', 'iupac', 'mw', 'name'},
            DataSource.CACTUS: {'smiles', 'cas', 'iupac', 'mw', 'name'},
            DataSource.CHEMBL: {'smiles', 'iupac', 'mw', 'name'},
            DataSource.CAS: {'smiles', 'cas', 'mw', 'name'}
        }

        self._supported_id_types = {
            DataSource.PUBCHEM: {'smiles', 'cas', 'name'},
            DataSource.CHEMSPIDER: {'smiles', 'cas', 'name'},
            DataSource.COMPTOX: {'smiles', 'cas', 'name'},
            DataSource.CACTUS: {'smiles', 'cas', 'name'},
            DataSource.CHEMBL: {'smiles', 'name'},
            DataSource.CAS: {'smiles', 'cas', 'name'}
        }
    def _increment_request_count(self):
        """
        Increment the request count and check if a break is needed.
        If the number of requests exceeds the batch limit, pause the execution for a specified duration.
        """
        if self.batch_limit <= 0:
            return
        with self._counter_lock:
            self.request_count += 1
            if self.request_count >= self.batch_limit:
                logger.info(f"Limit reached. Pausing for {self.rest_duration}s...")
                for remaining in range(self.rest_duration, 0, -1):
                    if self.status_callback:
                        try:
                            msg = f"API Rate Limit Reached. Cooling down: **{remaining}s** remaining..."
                            self.status_callback(msg)
                        except:
                            pass
                    time.sleep(1)
                self.request_count = 0
                if self.status_callback:
                    try:
                        self.status_callback("Cooldown finished. Resuming...")
                    except:
                        pass

    def get_properties_batch(self,
                             identifiers: List[str],
                             properties: Set[str],
                             identifier_type: str,
                             progress_callback: Optional[Callable[[int, int], None]] = None) \
            -> List[Dict[str, Optional[str]]]:
        """Batch acquisition of attributes, using multiple threads internally."""

        def is_valid_id(identifier):
            if identifier is None:
                return False
            if pd.isna(identifier):
                return False
            if isinstance(identifier, str) and not identifier.strip():
                return False
            return True

        unique_identifiers = list(set([i for i in identifiers if is_valid_id(i)]))
        results_map = {}

        total_tasks = len(unique_identifiers)
        completed_count = 0

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_identifier = {
                executor.submit(self.get_single_property, i, properties, identifier_type): i
                for i in unique_identifiers
            }

            for future in tqdm(as_completed(future_to_identifier), total=len(future_to_identifier),
                               desc=f"Query by {identifier_type}"):
                identifier = future_to_identifier[future]
                try:
                    results_map[identifier] = future.result()
                except Exception as e:
                    logger.error(f"An error occurred while processing the identifier '{identifier}': {e}")
                    error_res: Dict[str, Optional[str]] = {prop: None for prop in properties}
                    error_res['Data_Source'] = f"Thread Crash: {type(e).__name__}"
                    results_map[identifier] = error_res
                completed_count += 1
                if progress_callback:
                    progress_callback(completed_count, total_tasks)

        final_results = []
        for original_id in identifiers:
            if is_valid_id(original_id) and original_id in results_map:
                final_results.append(results_map[original_id])
            else:
                empty_res: Dict[str, Optional[str]] = {prop: None for prop in properties}
                if not is_valid_id(original_id):
                    empty_res['Data_Source'] = 'Invalid Identifier'
                else:
                    empty_res['Data_Source'] = 'Processing Failed (Missing from Results)'
                final_results.append(empty_res)
        return final_results

    def get_single_property(self,
                            identifier: str,
                            properties: Set[str],
                            identifier_type: str) -> Dict[str, Optional[str]]:
        """Query a single identifier in order of priority."""
        final_result: Dict[str, Optional[str]] = {prop: None for prop in properties}
        final_result['Data_Source'] = None
        needed_props = properties.copy()
        source_to_fields: Dict[str, List[str]] = {}
        source_errors = []

        for source in self.sources:
            if not needed_props:
                break
            fetch_fn = self._fetch_functions.get(source)
            if not fetch_fn:
                continue

            if identifier_type not in self._supported_id_types.get(source, set()):
                source_errors.append(f"{source}: Unsupported Input Type ('{identifier_type}')")
                continue

            supported_for_source = self._supported_props.get(source, set())
            props_to_request = needed_props.intersection(supported_for_source)

            if not props_to_request:
                source_errors.append(f"{source}: Unsupported Output Props")
                continue

            try:
                source_result = self._retry_wrapper(fetch_fn, identifier, needed_props, identifier_type)

                if not source_result:
                    source_errors.append(f"{source}: Not Found")
                    continue

                found_this_round = set()
                for prop, value in source_result.items():
                    if prop in needed_props and value is not None:
                        final_result[prop] = value
                        found_this_round.add(prop)

                if found_this_round:
                    source_to_fields[source] = list(found_this_round)
                    needed_props -= found_this_round

                missing_but_supported = props_to_request - found_this_round
                if missing_but_supported:
                    source_errors.append(f"{source}: Missing {list(missing_but_supported)}")

            except Exception as e:
                err_msg = str(e)
                source_errors.append(f"{source}: {err_msg}")
                logger.debug(f"Data source '{source}' failed to query '{identifier}': {e}. Trying the next one.")

        if source_to_fields:
            parts = [f"{src} ({', '.join(sorted(fields))})" for src, fields in source_to_fields.items()]
            final_result['Data_Source'] = " | ".join(parts)
        else:
            if source_errors:
                final_result['Data_Source'] = "Failed -> " + " | ".join(source_errors)
            else:
                final_result['Data_Source'] = "Failed -> Unsupported Query Strategy"

        return final_result

    def _retry_wrapper(self,
                       func: Callable,
                       identifier: str,
                       properties: Set[str],
                       identifier_type: str) -> Dict[str, Any]:
        """Wrapper that implements retry logic for failed requests."""
        last_error = None
        source_name = func.__name__.replace('_fetch_via_', '')
        for attempt in range(self.retries):
            try:
                time.sleep(self.interval)
                return func(identifier, properties, identifier_type)

            except Exception as e:
                last_error = str(e) or type(e).__name__
                error_msg_lower = last_error.lower()

                is_not_found = ((isinstance(e, requests.exceptions.HTTPError) and e.response.status_code in [400,404])
                                or any(kw in error_msg_lower for kw in
                                       ['400', '404', 'not found', 'invalid query', 'no record']))

                if is_not_found:
                    return {}

                is_auth_error = ((isinstance(e, requests.exceptions.HTTPError) and e.response.status_code in [401,403])
                                 or any(kw in error_msg_lower for kw in
                                        ['api key', 'unauthorized', 'forbidden', '401', '403']))

                if is_auth_error:
                    raise RuntimeError(f"Auth Error: {last_error}")

                if '429' in error_msg_lower or 'rate limit' in error_msg_lower:
                    last_error = "Rate Limit Exceeded (429)"
                elif any(kw in error_msg_lower for kw in ['500', '502', '503', '504', 'server error']):
                    import re
                    match = re.search(r'(50[0-4])', last_error)
                    code = match.group(1) if match else "5xx"
                    last_error = f"Server Error ({code})"
                elif any(kw in error_msg_lower for kw in ['timeout', 'timed out']):
                    last_error = "Timeout"
                elif 'connection' in error_msg_lower:
                    last_error = "Connection Error"

            if attempt < self.retries - 1:
                sleep_time = self.delay * (attempt + 1)
                logger.warning(
                    f"[{source_name}] Query '{identifier[:30]}...' failed ({last_error}). "
                    f"Retrying {attempt + 1}/{self.retries - 1} in {sleep_time}s..."
                )
                time.sleep(sleep_time)
            else:
                logger.debug(f"All {self.retries} attempts failed for '{identifier}'")

        raise RuntimeError(last_error)

    def _validate_and_clean_cas(self, cas: str) -> Optional[str]:
        if pd.isna(cas):
            return None

        cas_str = str(cas).strip()
        if not cas_str or cas_str.lower() in ['nan', 'none', 'null', 'n/a']:
            return None

        cas_str = re.sub(r'^(casrn|cas#|cas)[\s:]*', '', cas_str, flags=re.IGNORECASE).strip()
        if re.search(r'[a-zA-Z]', cas_str):
            return None
        cas_str = re.sub(r'[^\d\-]', '', cas_str)
        if '-' not in cas_str and len(cas_str) >= 5:
            cas_str = f"{cas_str[:-3]}-{cas_str[-3:-1]}-{cas_str[-1]}"

        if re.match(r'^\d{2,7}-\d{2}-\d$', cas_str) and self._validate_cas_checksum(cas_str):
            return cas_str
        return None

    @staticmethod
    def _validate_cas_checksum(cas: str) -> bool:
        """Verify the check digit of the CAS number."""
        try:
            digits = cas.replace('-', '')
            check_digit = int(digits[-1])
            s = sum(int(digit) * (i + 1) for i, digit in enumerate(digits[-2::-1]))
            return s % 10 == check_digit
        except (ValueError, IndexError):
            return False

    def _fetch_via_pubchem_sdk(self,
                               identifier: str,
                               properties: Set[str],
                               id_type: str) -> Dict[str, Optional[str]]:
        id_type_map = {'cas': 'name', 'name': 'name', 'smiles': 'smiles'}
        if id_type not in id_type_map:
            return {}
        self._increment_request_count()
        current_version = pubchempy.__version__
        result = {}

        if parse_version(current_version) < parse_version("1.0.5"):
            if not hasattr(self, '_pubchem_version_warned'):
                logger.warning(
                    f"Detected PubChemPy version {current_version}. "
                    f"It is strongly recommended to update to v1.0.5+ (`pip install pubchempy --upgrade`)"
                )
                self._pubchem_version_warned = True

            try:
                compounds = pubchempy.get_compounds(identifier, id_type_map[id_type])
            except Exception as e:
                raise RuntimeError(f"PubChem SDK error: {str(e)}")

            result = {}
            if compounds:
                c = compounds[0]
                if 'smiles' in properties:
                    result['smiles'] = getattr(c, 'canonical_smiles', None) or getattr(c, 'isomeric_smiles', None)
                if 'iupac' in properties:
                    result['iupac'] = getattr(c, 'iupac_name', None) or getattr(c, 'preferred_iupac_name', None)
                if 'mw' in properties:
                    result['mw'] = c.molecular_weight
                if c.synonyms is not None:
                    if 'cas' in properties:
                        result['cas'] = next((syn for syn in c.synonyms if self._validate_and_clean_cas(syn)), None)
                    if 'name' in properties and c.synonyms:
                        result['name'] = c.synonyms[0]
            return result

        else:
            try:
                compounds = pubchempy.get_compounds(identifier, id_type_map[id_type])
            except Exception as e:
                raise RuntimeError(f"PubChem SDK error: {str(e)}")

            if compounds:
                c = compounds[0]
                if 'smiles' in properties:
                    val = getattr(c, 'smiles', None)
                    if not val:
                        val = getattr(c, 'connectivity_smiles', None)
                    result['smiles'] = val

                if 'iupac' in properties:
                    result['iupac'] = getattr(c, 'iupac_name', None)

                if 'mw' in properties:
                    result['mw'] = getattr(c, 'molecular_weight', None)

                if 'cas' in properties or 'name' in properties:
                    syns = c.synonyms
                    if syns:
                        if 'name' in properties:
                            result['name'] = syns[0]
                        if 'cas' in properties:
                            result['cas'] = next((syn for syn in syns if self._validate_and_clean_cas(syn)), None)

            return result

    def _fetch_via_pubchem_api(self,
                               identifier: str,
                               properties: Set[str],
                               id_type: str) -> Dict[str, Optional[str]]:
        id_type_map = {'cas': 'name', 'name': 'name', 'smiles': 'smiles'}
        if id_type not in id_type_map:
            return {}
        base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        single_value_props = {'CanonicalSMILES', 'ConnectivitySMILES', 'IUPACName', 'MolecularWeight'}
        prop_str = ",".join(single_value_props)
        final_result = {}

        def get_and_parse(id_value: str, id_namespace: str) -> Optional[Dict]:
            parsed_result = {}

            prop_url = f"{base_url}/compound/{id_namespace}/{quote(id_value)}/property/{prop_str}/JSON"
            self._increment_request_count()
            response = self.session.get(prop_url, timeout=15)

            if response.status_code == 404:
                return None
            response.raise_for_status()

            data = response.json()
            if props := data.get('PropertyTable', {}).get('Properties', []):
                vals = props[0]
                if 'smiles' in properties: parsed_result['smiles'] = vals.get('CanonicalSMILES') or vals.get(
                    'ConnectivitySMILES')
                if 'iupac' in properties: parsed_result['iupac'] = vals.get('IUPACName')
                if 'mw' in properties: parsed_result['mw'] = vals.get('MolecularWeight')

            if 'name' in properties or 'cas' in properties:
                syn_url = f"{base_url}/compound/{id_namespace}/{quote(id_value)}/synonyms/JSON"
                self._increment_request_count()
                syn_response = self.session.get(syn_url, timeout=15)
                if syn_response.status_code != 404:
                    syn_response.raise_for_status()
                    syn_data = syn_response.json()
                    if synonyms := syn_data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym'):
                        if 'name' in properties:
                            parsed_result['name'] = synonyms[0]
                        if 'cas' in properties:
                            parsed_result['cas'] = next((syn for syn in synonyms if self._validate_and_clean_cas(syn)),
                                                        None)

            return parsed_result if parsed_result else None

        result = get_and_parse(identifier, id_type_map[id_type])
        if result:
            return result

        cid_url = f"{base_url}/compound/{id_type_map[id_type]}/{quote(identifier)}/cids/JSON"
        self._increment_request_count()
        cid_response = self.session.get(cid_url, timeout=15)
        if cid_response.status_code == 404:
            return {}
        cid_response.raise_for_status()

        if cids := cid_response.json().get('IdentifierList', {}).get('CID'):
            result = get_and_parse(str(cids[0]), 'cid')
            if result:
                return result

        return {}

    def _fetch_via_chemspider_sdk(self,
                                  identifier: str,
                                  properties: Set[str],
                                  id_type: str) -> Dict[str, Optional[str]]:
        if not self._chemspider_api_key:
            raise ValueError("ChemSpider requires an API key.")

        self._increment_request_count()
        cs = ChemSpider(self._chemspider_api_key)

        if id_type in ['cas', 'name']:
            query_id = cs.filter_name(identifier)
        elif id_type == 'smiles':
            query_id = cs.filter_smiles(identifier)
        else:
            return {}

        if not query_id:
            return {}

        for _ in range(3):
            status_info = cs.filter_status(query_id)
            status = status_info.get('status')
            if status == 'Complete':
                break
            if status == 'Failed':
                raise RuntimeError("ChemSpider SDK Filter Task Failed")
            time.sleep(2)
        else:
            raise RuntimeError("ChemSpider SDK Task Timeout")

        results_csids = cs.filter_results(query_id)
        if not results_csids:
            return {}
        csid = results_csids[0]

        field_map_request = {'smiles': 'SMILES', 'mw': 'MolecularWeight', 'name': 'CommonName'}
        fields_to_get = [field_map_request[prop] for prop in properties if prop in field_map_request]
        if not fields_to_get:
            return {}

        detail_data = cs.get_details(csid, fields=fields_to_get)

        result = {}
        if 'smiles' in properties:
            result['smiles'] = detail_data.get('smiles')
        if 'mw' in properties:
            result['mw'] = detail_data.get('molecularWeight')
        if 'name' in properties:
            result['name'] = detail_data.get('commonName')
        return result

    def _fetch_via_chemspider_api(self,
                                  identifier: str,
                                  properties: Set[str],
                                  id_type: str) -> Dict[str, Optional[str]]:
        if not self._chemspider_api_key:
            raise ValueError("ChemSpider requires an API key.")

        if id_type in ['cas', 'name']:
            filter_url = "https://api.rsc.org/compounds/v1/filter/name"
            payload = {"name": identifier}
        elif id_type == 'smiles':
            filter_url = "https://api.rsc.org/compounds/v1/filter/smiles"
            payload = {"smiles": identifier}
        else:
            return {}

        headers_post = {
            "apikey": self._chemspider_api_key,
            "Content-Type": "application/json",
            "Accept": "application/json"
        }

        self._increment_request_count()
        response_post = self.session.post(filter_url, json=payload, headers=headers_post, timeout=15)
        response_post.raise_for_status()

        query_id = response_post.json().get('queryId')
        if not query_id:
            return {}

        status_url = f"https://api.rsc.org/compounds/v1/filter/{query_id}/status"
        headers_get = {
            "apikey": self._chemspider_api_key,
            "Accept": "application/json"
        }

        # 轮询状态
        for _ in range(3):
            self._increment_request_count()
            status_response = self.session.get(status_url, headers=headers_get, timeout=15)
            status_response.raise_for_status()
            status = status_response.json().get('status')

            if status == 'Complete':
                break
            if status == 'Failed':
                raise RuntimeError("ChemSpider Filter Task Failed")
            time.sleep(2)
        else:
            raise RuntimeError("ChemSpider Task Timeout")

        results_url = f"https://api.rsc.org/compounds/v1/filter/{query_id}/results"
        self._increment_request_count()
        response_get_results = self.session.get(results_url, headers=headers_get, timeout=15)
        response_get_results.raise_for_status()

        results_data = response_get_results.json()
        if not results_data.get('results'):
            return {}
        csid = results_data['results'][0]

        field_map = {'smiles': 'SMILES', 'mw': 'MolecularWeight', 'name': 'CommonName'}
        fields_to_get = ",".join([field_map[prop] for prop in properties if prop in field_map])
        if not fields_to_get:
            return {}

        detail_url = f"https://api.rsc.org/compounds/v1/records/{csid}/details?fields={fields_to_get}"

        self._increment_request_count()
        response_get_details = self.session.get(detail_url, headers=headers_get, timeout=15)
        response_get_details.raise_for_status()
        detail_data = response_get_details.json()

        result = {}
        if 'smiles' in properties:
            result['smiles'] = detail_data.get('smiles')
        if 'mw' in properties:
            result['mw'] = detail_data.get('molecularWeight')
        if 'name' in properties:
            result['name'] = detail_data.get('commonName')
        return result

    @staticmethod
    def _get_inchikey_from_smiles(smiles: str) -> Optional[str]:
        if not smiles or pd.isna(smiles):
            return None
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                logger.debug(f"Unable to parse SMILES to generate InChIKey: {smiles}")
                return None
            return Chem.MolToInchiKey(mol)
        except Exception as e:
            return None

    @staticmethod
    def _get_inchi_from_smiles(smiles: str) -> Optional[str]:
        if not smiles or pd.isna(smiles):
            return None
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                logger.debug(f"Unable to parse SMILES to generate InChI: {smiles}")
                return None
            return Chem.MolToInchi(mol)
        except Exception as e:
            return None
        finally:
            lg.setLevel(RDLogger.INFO)

    def _fetch_via_comptox_sdk(self,
                               identifier: str,
                               properties: Set[str],
                               id_type: str) -> Dict[str, Optional[str]]:
        if not self._comptox_api_key:
            raise ValueError("CompTox requires an API key.")

        if identifier is None:
            return {}

        search_identifier = identifier
        if id_type == 'smiles':
            search_identifier = self._get_inchikey_from_smiles(identifier)
            if not search_identifier:
                return {}

        comptox = ctx.Chemical(x_api_key=self._comptox_api_key)

        self._increment_request_count()
        search_results = comptox.search(by='equals', query=search_identifier)

        if not search_results or not isinstance(search_results, list):
            return {}

        dtxsid = search_results[0].get('dtxsid')
        if not dtxsid:
            return {}

        self._increment_request_count()
        detail_resp = comptox.details(by='dtxsid', query=dtxsid)
        result = {}
        if isinstance(detail_resp, list):
            if len(detail_resp) > 0:
                detail = detail_resp[0]
            else:
                return {}
        elif isinstance(detail_resp, dict):
            detail = detail_resp
        else:
            return {}

        if 'smiles' in properties:
            result['smiles'] = detail.get('smiles')
        if 'iupac' in properties:
            result['iupac'] = detail.get('iupacName')
        if 'mw' in properties:
            result['mw'] = detail.get('averageMass')
        if 'cas' in properties:
            result['cas'] = detail.get('casrn')
        if 'name' in properties:
            result['name'] = detail.get('preferredName')

        return result

    def _fetch_via_comptox_api(self,
                               identifier: str,
                               properties: Set[str],
                               id_type: str) -> Dict[str, Optional[str]]:
        if not self._comptox_api_key:
            raise ValueError("CompTox requires an API key.")

        search_base_url = "https://comptox.epa.gov/ctx-api/chemical/search/equal"
        search_identifier = identifier

        if id_type == 'smiles':
            search_identifier = self._get_inchikey_from_smiles(identifier)
            if not search_identifier:
                return {}
        elif id_type not in ['cas', 'name']:
            return {}
        search_url = f"{search_base_url}/{quote(search_identifier)}"
        headers = {
            'x-api-key': self._comptox_api_key
        }
        self._increment_request_count()
        response = self.session.get(search_url, headers=headers, timeout=15)

        if response.status_code == 404:
            return {}
        response.raise_for_status()

        search_results = response.json()
        if not search_results or not isinstance(search_results, list):
            return {}

        dtxsid = search_results[0].get('dtxsid')
        if not dtxsid:
            logger.debug(f"CompTox Search for '{identifier}' succeeded but returned no DTXSID.")
            return {}

        details_base_url = "https://comptox.epa.gov/ctx-api/chemical/detail/search/by-dtxsid"
        details_url = f"{details_base_url}/{quote(dtxsid)}"

        self._increment_request_count()
        details_response = self.session.get(details_url, headers=headers, timeout=15)
        if details_response.status_code == 404:
            return {}
        details_response.raise_for_status()

        chemical_data = details_response.json()
        if not chemical_data or not isinstance(chemical_data, dict):
            return {}

        result = {}
        if 'smiles' in properties:
            result['smiles'] = chemical_data.get('smiles')
        if 'iupac' in properties:
            result['iupac'] = chemical_data.get('iupacName')
        if 'mw' in properties:
            result['mw'] = chemical_data.get('averageMass')
        if 'cas' in properties:
            result['cas'] = chemical_data.get('casrn')
        if 'name' in properties:
            result['name'] = chemical_data.get('preferredName')
        return result

    def _fetch_via_cactus_api(self,
                              identifier: str,
                              properties: Set[str],
                              id_type: str) -> Dict[str, Optional[str]]:
        search_identifier = identifier
        if id_type == 'smiles':
            search_identifier = self._get_inchikey_from_smiles(identifier)
            if not search_identifier:
                return {}
        elif id_type not in ['cas', 'name']:
            return {}
        base_url = f"https://cactus.nci.nih.gov/chemical/structure/{quote(search_identifier)}"
        prop_map = {'cas': 'cas', 'iupac': 'iupac_name', 'smiles': 'smiles', 'mw': 'mw', 'name': 'names'}
        result = {}
        for prop in properties:
            if url_suffix := prop_map.get(prop):
                self._increment_request_count()
                resp = self.session.get(f"{base_url}/{url_suffix}", timeout=10)
                if resp.status_code == 404:
                    continue
                resp.raise_for_status()

                if resp.ok and resp.text:
                    text_content = resp.text.strip()
                    if not text_content:
                        continue
                    if prop == 'name':
                        result[prop] = text_content.split('\n')[0]
                    elif prop == 'cas':
                        potential_cas_list = text_content.split('\n')
                        valid_cas_list = [
                            cleaned_cas for cas in potential_cas_list
                            if (cleaned_cas := self._validate_and_clean_cas(cas))
                        ]
                        if valid_cas_list:
                            result[prop] = min(valid_cas_list, key=len)
                    else:
                        result[prop] = text_content
        return result

    def _fetch_via_chembl_api(self,
                              identifier: str,
                              properties: Set[str],
                              id_type: str) -> Dict[str, Optional[str]]:
        if id_type == 'cas':
            return {}
        base_url = "https://www.ebi.ac.uk/chembl/api/data/molecule"
        search_url = ""
        if id_type == 'smiles':
            search_url = f"{base_url}.json?molecule_structures__canonical_smiles__exact={quote(identifier)}"
        elif id_type == 'name':
            search_url = f"{base_url}.json?pref_name__iexact={quote(identifier)}"
        if not search_url: return {}
        self._increment_request_count()
        response = self.session.get(search_url, timeout=15)
        response.raise_for_status()
        data = response.json()
        result = {}
        if mols := data.get('molecules'):
            mol = mols[0]
            props = mol.get('molecule_properties') or {}
            structs = mol.get('molecule_structures') or {}
            if 'smiles' in properties:
                result['smiles'] = structs.get('canonical_smiles')
            if 'iupac' in properties and mol.get('pref_name'):
                result['iupac'] = mol.get('pref_name')
            if 'mw' in properties:
                result['mw'] = props.get('full_mwt')
            if 'name' in properties and mol.get('pref_name'):
                result['name'] = mol.get('pref_name')
        return result

    def _fetch_via_cas_api(self,
                           identifier: str,
                           properties: Set[str],
                           id_type: str) -> Dict[str, Optional[str]]:
        if not self._cas_api_key:
            raise ValueError("CAS Common Chemistry API requires an API key.")
        if id_type not in ['smiles', 'cas', 'name']:
            return {}
        base_url = "https://commonchemistry.cas.org/api"
        headers = {
            'X-Api-Key': self._cas_api_key
        }
        search_url = f"{base_url}/search"
        params = {'q': identifier}

        self._increment_request_count()

        search_response = self.session.get(search_url, params=params, headers=headers, timeout=15)
        if search_response.status_code == 404:
            return {}
        search_response.raise_for_status()

        search_data = search_response.json()
        if not search_data.get('count', 0) > 0 or not search_data.get('results'):
            return {}

        first_result = search_data['results'][0]
        cas_rn = first_result.get('rn')
        if not cas_rn:
            return {}

        detail_url = f"{base_url}/detail"
        detail_params = {'cas_rn': cas_rn}

        self._increment_request_count()
        detail_response = self.session.get(detail_url, params=detail_params, headers=headers, timeout=15)

        if detail_response.status_code == 404:
            return {}
        detail_response.raise_for_status()

        detail_data = detail_response.json()
        if not detail_data:
            return {}

        result = {}
        if 'smiles' in properties:
            result['smiles'] = detail_data.get('canonicalSmile') or detail_data.get('smile')

        if 'mw' in properties:
            result['mw'] = detail_data.get('molecularMass')

        if 'cas' in properties:
            result['cas'] = detail_data.get('rn')

        if 'name' in properties:
            result['name'] = detail_data.get('name')
            if not result['name'] and detail_data.get('synonyms'):
                result['name'] = detail_data['synonyms'][0]

        return result