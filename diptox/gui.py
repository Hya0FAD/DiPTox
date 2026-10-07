"""NiceGUI web interface for DiPTox.

The interface is a long-lived single-page application. Navigation never reruns
the application, and pipeline mutations are performed on a copy in a worker
thread before being committed atomically.
"""

from __future__ import annotations

import asyncio
import contextvars
import copy
import functools
import json
import multiprocessing
import os
import re
import shutil
import tempfile
import threading
import traceback
import uuid
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

os.environ["DIPTOX_GUI_MODE"] = "true"

import pandas as pd
from nicegui import events, ui

from .core import DiptoxPipeline
from .unit_processor import UnitProcessor
from . import user_reg


NONE_TEXT = "Not mapped"
FULL_NAME = "Data Integration and Processing for Computational Toxicology"
BRAND_IMAGE = Path(__file__).with_name("logo.svg")
PAGE_DEFINITIONS = (
    ("load", "Data loading", "upload_file"),
    ("preprocess", "Preprocessing", "science"),
    ("web", "Web requests", "language"),
    ("units", "Unit standardization", "straighten"),
    ("columns", "Column adjustments", "tune"),
    ("deduplicate", "Deduplication", "difference"),
    ("search", "Search & filter", "manage_search"),
    ("export", "Export", "download"),
)
SUPPORTED_SUFFIXES = {".csv", ".xlsx", ".xls", ".txt", ".sdf", ".mol", ".smi"}
HEADER_OPTION_SUFFIXES = {".csv", ".txt", ".smi"}


async def _to_thread(function: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    """Run blocking work with its context on Python 3.8 as well as newer Python."""
    context = contextvars.copy_context()
    operation = functools.partial(context.run, function, *args, **kwargs)
    return await asyncio.get_running_loop().run_in_executor(None, operation)


def _optional_text(value: Any) -> Optional[str]:
    text = "" if value is None else str(value).strip()
    return None if not text or text == NONE_TEXT else text


def _clone_pipeline(source: DiptoxPipeline) -> DiptoxPipeline:
    """Copy mutable pipeline state without copying a WebService thread lock."""
    clone = copy.copy(source)
    clone.df = source.df.copy(deep=True) if source.df is not None else None
    clone.excluded_df = source.excluded_df.copy(deep=True)
    clone.chem_processor = copy.deepcopy(source.chem_processor)
    clone.data_handler = copy.deepcopy(source.data_handler)
    clone.deduplicator = copy.deepcopy(source.deduplicator)
    clone.web_service = None
    clone._audit_log = copy.deepcopy(source._audit_log)
    clone._history = []
    for checkpoint in source._history:
        item = copy.deepcopy(
            {key: value for key, value in checkpoint.items() if key not in {"df", "excluded_df"}}
        )
        item["df"] = checkpoint["df"].copy(deep=True)
        item["excluded_df"] = checkpoint.get("excluded_df", pd.DataFrame()).copy(deep=True)
        clone._history.append(item)
    clone._dedup_unit_settings = copy.deepcopy(source._dedup_unit_settings)
    clone._current_dedup_config = copy.deepcopy(source._current_dedup_config)
    clone._structure_derived_columns = set(source._structure_derived_columns)
    clone._source_index_labels = copy.deepcopy(source._source_index_labels)
    clone._row_lineage = copy.deepcopy(source._row_lineage)
    clone.source_columns = list(source.source_columns)
    clone._input_column_aliases = dict(source._input_column_aliases)
    return clone


class BusyError(RuntimeError):
    """Raised when a second operation is submitted while one is running."""


class PipelineStore:
    """Thread-safe job state and atomic ownership of a DiPTox pipeline."""

    def __init__(self, pipeline: Optional[DiptoxPipeline] = None):
        self.pipeline = pipeline or DiptoxPipeline()
        self.busy = False
        self.description = "Ready"
        self.status = "No operation is running"
        self.current = 0
        self.total = 0
        self.last_error: Optional[str] = None
        self._state_lock = threading.Lock()

    def snapshot(self) -> tuple[bool, str, str, int, int, Optional[str]]:
        with self._state_lock:
            return self.busy, self.description, self.status, self.current, self.total, self.last_error

    def _begin(self, description: str) -> None:
        with self._state_lock:
            if self.busy:
                raise BusyError("Another operation is still running. You can keep browsing while it finishes.")
            self.busy = True
            self.description = description
            self.status = description
            self.current = 0
            self.total = 0
            self.last_error = None

    def _set_progress(self, current: int, total: int) -> None:
        with self._state_lock:
            self.current = max(0, int(current))
            self.total = max(0, int(total))

    def _set_status(self, status: str) -> None:
        with self._state_lock:
            self.status = str(status)

    def _finish(self, error: Optional[str] = None) -> None:
        with self._state_lock:
            self.busy = False
            self.last_error = error
            if error:
                self.status = "Operation failed"
            else:
                self.description = "Ready"
                self.status = "No operation is running"
                self.current = 0
                self.total = 0

    async def run_pipeline_job(
        self,
        description: str,
        action: Callable[[DiptoxPipeline, Callable[[int, int], None], Callable[[str], None]], None],
    ) -> DiptoxPipeline:
        """Run a mutation on a clone and commit only after successful completion."""
        self._begin(description)
        source = self.pipeline

        def operation() -> DiptoxPipeline:
            working = _clone_pipeline(source)
            action(working, self._set_progress, self._set_status)
            return working

        try:
            result = await _to_thread(operation)
            self.pipeline = result
            self._finish()
            return result
        except Exception as exc:
            self._finish(traceback.format_exc())
            raise exc

    async def run_read_job(
        self,
        description: str,
        operation: Callable[[Callable[[int, int], None], Callable[[str], None]], Any],
    ) -> Any:
        """Run a non-mutating operation without blocking the NiceGUI event loop."""
        self._begin(description)
        try:
            result = await _to_thread(operation, self._set_progress, self._set_status)
            self._finish()
            return result
        except Exception as exc:
            self._finish(traceback.format_exc())
            raise exc


@dataclass
class InterfaceState:
    upload_path: Optional[Path] = None
    upload_name: str = ""
    detected_columns: list[str] = field(default_factory=list)
    sheet_name: Optional[str] = None
    unit_rules: dict[tuple[str, str], str] = field(default_factory=dict)


def _safe_cell(value: Any) -> Any:
    if value is None or value is pd.NA:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        try:
            value = value.item()
        except (ValueError, TypeError):
            pass
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if hasattr(value, "isoformat"):
        try:
            return value.isoformat()
        except (TypeError, ValueError):
            pass
    return str(value)


def _grid_options(frame: pd.DataFrame, limit: int = 200) -> dict[str, Any]:
    """Build JSON-safe AG Grid options from a DataFrame preview."""
    subset = frame.head(limit)
    column_defs = [
        {"field": f"c{index}", "headerName": str(column), "minWidth": 130}
        for index, column in enumerate(subset.columns)
    ]
    rows = [
        {f"c{index}": _safe_cell(value) for index, value in enumerate(row)}
        for row in subset.itertuples(index=False, name=None)
    ]
    return {
        "columnDefs": column_defs,
        "rowData": rows,
        "defaultColDef": {"sortable": True, "filter": True, "resizable": True},
        "pagination": True,
        "paginationPageSize": 25,
        "animateRows": False,
    }


@contextmanager
def _work_card(title: str, description: str = ""):
    with ui.card().classes("work-card w-full") as card:
        ui.label(title).classes("section-title")
        if description:
            ui.label(description).classes("section-copy")
        yield card


class DiptoxWebApp:
    """One browser client's persistent DiPTox workbench."""

    MAPPINGS = (
        ("SMILES column", "smiles"),
        ("Name column", "name"),
        ("CAS column", "cas"),
        ("ID column", "id"),
        ("Target column", "target"),
        ("Unit column", "unit"),
    )

    def __init__(self):
        self.store = PipelineStore()
        self.state = InterfaceState()
        self.mapping: dict[str, Any] = {}
        self.unit_rule_inputs: dict[tuple[str, str], Any] = {}
        self._export_selection_mode = 'recommended'
        self._export_refreshing = False
        self._column_value_options = []
        self._build()
        self._refresh_data_controls()
        self._build_registration_dialog()

    def _build_registration_dialog(self) -> None:
        if user_reg.is_registered_or_skipped():
            return
        with ui.dialog().props("persistent") as dialog, ui.card().classes("work-card w-[520px] max-w-[94vw]"):
            with ui.row().classes("w-full items-center no-wrap"):
                ui.label("Welcome to DiPTox").classes("text-xl font-bold")
                ui.space()
                ui.button(icon="close", on_click=dialog.close).props(
                    'flat round dense aria-label="Close questionnaire for now"'
                ).tooltip("Close for now. Ask again next time.")
            ui.label("Optional user questionnaire").classes("section-title")
            ui.label(
                "Help the DiPTox developers understand who uses the tool. Submitting sends only the name, "
                "affiliation and optional email entered below via Google Forms or the Feishu backup channel. "
                "Your molecular datasets are not included. You can skip without affecting any features."
            ).classes("section-copy")
            name = ui.input("Name").props("outlined dense").classes("w-full")
            affiliation = ui.input("Affiliation / unit").props("outlined dense").classes("w-full")
            email = ui.input("Email (optional)").props("outlined dense type=email").classes("w-full")
            feedback = ui.label("").classes("text-sm text-slate-600")

            async def skip() -> None:
                await _to_thread(user_reg.save_status, "skipped")
                if not user_reg.is_registered_or_skipped():
                    ui.notify("Could not save your choice. The questionnaire may appear next time.", type="warning")
                dialog.close()

            async def submit() -> None:
                values = tuple((entry.value or "").strip() for entry in (name, affiliation, email))
                if not values[0] or not values[1]:
                    feedback.set_text("Enter your name and affiliation, or choose Skip.")
                    return
                submit_button.disable()
                skip_button.disable()
                submit_button.props("loading")
                feedback.set_text("Submitting your questionnaire…")
                try:
                    success, message = await _to_thread(user_reg.submit_info, *values)
                    if success:
                        if not user_reg.is_registered_or_skipped():
                            ui.notify("Submitted, but the local preference could not be saved.", type="warning")
                        else:
                            ui.notify(message, type="positive")
                        dialog.close()
                    else:
                        feedback.set_text(f"{message} You can retry or skip.")
                except Exception:
                    feedback.set_text("Submission failed. You can retry or skip.")
                finally:
                    submit_button.props(remove="loading")
                    submit_button.enable()
                    skip_button.enable()

            with ui.row().classes("action-row"):
                skip_button = ui.button("Skip", on_click=skip).props("outline no-caps color=primary")
                submit_button = ui.button("Submit questionnaire", on_click=submit).props("unelevated no-caps color=primary").classes("primary-action")
        ui.timer(0.1, dialog.open, once=True)

    @property
    def pipeline(self) -> DiptoxPipeline:
        return self.store.pipeline

    def _build(self) -> None:
        ui.colors(primary="#0f766e", secondary="#334155", accent="#0d9488", positive="#15803d", negative="#b42318")
        ui.add_head_html(
            """
            <style>
            :root { --canvas:#f4f6f8; --surface:#fff; --ink:#172033; --muted:#657184;
              --line:#dbe1e8; --soft:#eef3f4; --accent:#0f766e; --accent-soft:#e8f4f2; }
            /* Quasar locks body for dialogs. Keep the scrollbar on that same element
               so opening a dialog does not introduce a second scrollbar gutter. */
            body { overflow-y:scroll; background:var(--canvas); color:var(--ink); font-family:Inter,"Segoe UI",sans-serif; }
            .q-layout,.q-page-container { background:var(--canvas); }
            .app-header { height:68px; background:rgba(255,255,255,.98); border-bottom:1px solid var(--line); }
            .brand-image { width:112px; height:48px; object-fit:contain; flex-shrink:0; }
            .brand-title { color:var(--ink); font-size:17px; line-height:1.05; font-weight:750; letter-spacing:-.02em; }
            .brand-subtitle { color:var(--muted); font-size:12px; line-height:1.35; max-width:470px; }
            .app-drawer { background:var(--surface); border-right:1px solid var(--line); }
            .nav-caption { color:#8a94a3; font-size:10px; font-weight:700; letter-spacing:.12em; text-transform:uppercase; }
            .vertical-nav .q-tabs__content { gap:4px; }
            .vertical-nav .q-tab { min-height:43px; padding:0 12px; border-radius:8px; justify-content:flex-start; color:#536174; }
            .vertical-nav .q-tab__content { flex-direction:row; justify-content:flex-start; gap:11px; min-width:100%; }
            .vertical-nav .q-tab__label { font-size:13px; font-weight:600; text-transform:none; letter-spacing:0; }
            .vertical-nav .q-tab--active { background:var(--accent-soft); color:var(--accent); }
            .vertical-nav .q-tab__indicator { display:none; }
            .page-shell { width:100%; max-width:1460px; margin:0 auto; padding:28px 32px 92px; }
            .page-panels,.page-panels .q-panel,.page-panels .q-tab-panel { background:transparent!important; }
            .page-panels .q-tab-panel { padding:0; }
            .page-title { color:var(--ink); font-size:28px; line-height:1.16; font-weight:760; letter-spacing:-.035em; }
            .work-card { background:var(--surface); border:1px solid var(--line); border-radius:12px;
              box-shadow:0 1px 2px rgba(20,31,48,.035); padding:20px; }
            .section-title { color:var(--ink); font-size:15px; font-weight:700; letter-spacing:-.012em; }
            .section-copy { color:var(--muted); font-size:12px; line-height:1.5; margin-top:-3px; margin-bottom:7px; }
            .grid-2 { display:grid!important; grid-template-columns:repeat(2,minmax(0,1fr)); gap:14px 18px; width:100%; }
            .grid-2 > .work-card { align-self:stretch; }
            .grid-3 { display:grid!important; grid-template-columns:repeat(3,minmax(0,1fr)); gap:14px 18px; width:100%; }
            .settings-grid { display:grid!important; grid-template-columns:repeat(3,minmax(150px,1fr)); gap:10px 16px; width:100%; }
            /* Shared form treatment, including controls inside dialogs. Keep the
               native Quasar labels and input associations for keyboard access. */
            .q-field--outlined { font-size:13px; }
            .q-field--outlined .q-field__control { border-radius:8px; background:#f8fafb;
              min-height:44px; color:var(--accent); }
            .q-field--outlined:not(.q-textarea) .q-field__control,
            .q-field--outlined .q-field__marginal { height:44px; }
            .q-field--outlined .q-field__control:before { border-color:#cbd5df; }
            .q-field--outlined:hover .q-field__control:before { border-color:#91aaa8; }
            .q-field--outlined.q-field--focused .q-field__control { background:#fff;
              box-shadow:0 0 0 3px rgba(15,118,110,.1); }
            .q-field--outlined.q-field--labeled { padding-top:23px; }
            .q-field--outlined .q-field__label { top:-22px; left:0; max-width:100%;
              transform:none!important; font-size:12px; line-height:18px; font-weight:600;
              color:#455468; transition:color .15s; }
            .q-field--outlined.q-field--focused .q-field__label { color:var(--accent); }
            .q-field--outlined .q-field__native,
            .q-field--outlined .q-field__input { padding-top:10px; padding-bottom:10px;
              min-height:44px; color:var(--ink); }
            .q-field--outlined .q-field__control-container { padding-top:0!important; }
            .q-field--outlined input::placeholder,
            .q-field--outlined textarea::placeholder { opacity:1; color:#718096; }
            .q-field--outlined.q-field--auto-height .q-field__control { height:auto; }
            .q-field--outlined.q-field--auto-height .q-field__native { padding:5px 0; }
            /* Search inputs sit inside the select's padded native container;
               they must not inherit the full-height standalone input sizing. */
            .q-field--outlined.q-select .q-field__input {
              min-height:24px; height:24px; padding-top:0; padding-bottom:0; }
            .q-field--outlined.q-textarea .q-field__native { padding:12px 0; line-height:1.65; }
            .q-field--outlined.q-field--error .q-field__control:before { border-color:#b42318; }
            .q-field--outlined.q-field--error .q-field__label { color:#b42318; }
            .q-field .q-chip { background:var(--accent-soft); color:#115e59;
              border:1px solid #cce5e1; border-radius:6px; font-size:12px; margin:3px 5px 3px 0; }
            .q-menu { border:1px solid var(--line); border-radius:10px;
              box-shadow:0 8px 24px rgba(20,31,48,.12); padding:5px; }
            .q-menu .q-item { min-height:38px; padding:9px 12px; border-radius:6px; font-size:13px; }
            .q-menu .q-item--active { background:var(--accent-soft); color:#115e59; font-weight:600; }
            .q-menu .q-item[aria-selected=true] .q-item__section--main:after {
              content:'✓'; position:absolute; right:12px; color:var(--accent); }
            .q-menu .q-item[aria-selected=true] { padding-right:34px; }
            .q-checkbox { min-height:36px; gap:4px; }
            .q-checkbox__inner { font-size:34px; color:#8999aa; }
            .q-checkbox__inner--truthy,.q-checkbox__inner--indet { color:var(--accent); }
            .q-checkbox__bg { border-radius:4px; border-width:1.5px; }
            .q-checkbox__label { font-size:13px; line-height:1.45; color:#354458; }
            .q-checkbox:focus-visible { outline:2px solid var(--accent); outline-offset:2px; border-radius:6px; }
            .q-btn { border-radius:8px; letter-spacing:0; font-weight:650; text-transform:none; }
            .primary-action { min-height:40px; padding:0 18px; }
            .action-row { width:100%; justify-content:flex-end; align-items:center; gap:12px;
              margin-top:auto; padding-top:16px; }
            .action-row .q-btn { min-height:40px; }
            .quiet-action { color:#425067; }
            .data-dialog .q-dialog__inner { padding:0; }
            .data-dialog .data-surface { width:100vw; max-width:100vw; height:68vh; max-height:85vh;
              overflow:auto; border-radius:16px 16px 0 0; border-bottom:0;
              box-shadow:0 -8px 32px rgba(20,31,48,.12); }
            .data-dialog .q-dialog__backdrop { background:transparent; }
            .metric { min-width:112px; padding:10px 12px; border-left:1px solid var(--line); }
            .metric-label { color:#8a94a3; font-size:10px; font-weight:700; text-transform:uppercase; letter-spacing:.08em; }
            .metric-value { color:var(--ink); font-size:14px; font-weight:700; }
            .empty-state { min-height:215px; display:flex!important; align-items:center; justify-content:center;
              flex-direction:column; color:var(--muted); border:1px dashed #cbd4df; border-radius:10px; background:#fafbfc; }
            .status-footer { min-height:50px; background:#fff; border-top:1px solid var(--line); padding:0 24px; }
            .status-footer { cursor:pointer; }
            .status-footer:hover { background:#f7faf9; }
            .status-footer:focus-visible { outline:2px solid var(--accent); outline-offset:-2px; }
            .status-dot { width:8px; height:8px; border-radius:50%; background:#9aa5b3; }
            .status-dot.busy { background:var(--accent); box-shadow:0 0 0 4px var(--accent-soft); }
            .upload-zone { width:100%; border-radius:10px; overflow:hidden; }
            .upload-zone .q-uploader__header { background:#f7faf9; color:var(--ink); border-bottom:1px solid var(--line); }
            .context-strip { background:#f7faf9; border:1px solid #d8e7e4; border-radius:9px; padding:11px 13px; color:#3d5958; }
            .rules-dialog-card { border:1px solid var(--line); border-radius:14px; }
            .rules-tabs { width:calc(100% - 48px); margin:0 24px; }
            .segmented-tabs { padding:4px;
              background:#eef2f4; border:1px solid var(--line); border-radius:10px; }
            .segmented-tabs .q-tabs__content { gap:4px; }
            .segmented-tabs .q-tab { min-height:38px; padding:0 14px; border-radius:7px; color:#536174; }
            .segmented-tabs .q-tab__label { font-size:13px; font-weight:600; text-transform:none; }
            .segmented-tabs .q-tab--active { background:#fff; color:var(--accent);
              box-shadow:0 1px 3px rgba(20,31,48,.1); }
            .segmented-tabs .q-tab__indicator { display:none; }
            .rules-panels .q-tab-panel { padding:0; }
            .rules-layout { display:grid!important; grid-template-columns:minmax(0,1.35fr) minmax(0,1fr);
              gap:24px; width:100%; align-items:start; }
            .rule-list { height:320px; overflow:auto; scrollbar-gutter:stable; width:100%;
              background:#f7f9fa; border:1px solid var(--line); border-radius:10px; padding:12px;
              display:grid!important; grid-template-columns:repeat(2,minmax(0,1fr)); gap:8px;
              align-content:start; }
            .rule-list--atoms { grid-template-columns:repeat(auto-fit,minmax(48px,1fr)); }
            .rule-list--neutral { grid-template-columns:minmax(0,1fr); }
            .rule-tile { min-width:0; min-height:40px; display:flex; align-items:center;
              padding:9px 11px; border:1px solid #dce5e9; border-radius:7px; background:#fff;
              color:#354458; font:12px/1.55 Consolas,monospace; overflow-wrap:anywhere; }
            .rule-list--atoms .rule-tile { justify-content:center; color:#115e59;
              background:var(--accent-soft); border-color:#cce5e1; font-size:14px; font-weight:600; }
            .rule-empty { grid-column:1/-1; color:var(--muted); padding:12px; font-size:13px; }
            @media (max-width:700px) { .rules-layout{grid-template-columns:minmax(0,1fr)}
              .rule-list{height:240px} .segmented-tabs .q-tab{padding:0 8px} }
            @media (max-width:980px) { .page-shell{padding:22px 18px 84px} .grid-2,.grid-3,.settings-grid{grid-template-columns:1fr} .header-metrics{display:none!important} }
            </style>
            """
        )

        with ui.header(fixed=True).classes("app-header items-center px-5 no-wrap"):
            with ui.row().classes("items-center gap-3 no-wrap"):
                ui.image(BRAND_IMAGE).props('fit=contain alt="DiPTox logo"').classes("brand-image")
                with ui.column().classes("gap-0"):
                    ui.label("DiPTox").classes("brand-title")
                    ui.label(FULL_NAME).classes("brand-subtitle")
            ui.space()
            with ui.row().classes("header-metrics items-center gap-0 no-wrap"):
                with ui.column().classes("metric gap-0"):
                    ui.label("Dataset").classes("metric-label")
                    self.header_dataset = ui.label("Not loaded").classes("metric-value")
                with ui.column().classes("metric gap-0"):
                    ui.label("Workspace").classes("metric-label")
                    self.header_stage = ui.label("Data loading").classes("metric-value")

        with ui.left_drawer(fixed=True, bordered=True).classes("app-drawer p-4"):
            self.nav_tabs = ui.tabs(on_change=self._navigation_changed).props("vertical align=left").classes("vertical-nav w-full")
            self.tabs: dict[str, Any] = {}
            with self.nav_tabs:
                for name, label, icon in PAGE_DEFINITIONS:
                    self.tabs[name] = ui.tab(name, label=label, icon=icon)
            ui.separator().classes("my-4")
            with ui.column().classes("w-full gap-2"):
                ui.label("Pipeline state").classes("nav-caption px-2")
                self.sidebar_summary = ui.label("No dataset loaded").classes("text-xs text-slate-500 px-2 leading-5")
                self.undo_button = ui.button("Undo last step", icon="undo", on_click=self._undo).props("flat no-caps").classes("quiet-action w-full")

        with ui.footer(fixed=True).classes("status-footer items-center no-wrap").props(
            'role=button tabindex=0 aria-label="Open workspace data" aria-haspopup=dialog'
        ).on("click", lambda: self.data_dialog.open()).on(
            "keydown.enter.prevent", lambda: self.data_dialog.open()
        ).on("keydown.space.prevent", lambda: self.data_dialog.open()):
            self.status_dot = ui.html('<div class="status-dot"></div>')
            with ui.column().classes("gap-0 ml-2 min-w-0"):
                self.job_title = ui.label("Ready").classes("text-xs font-semibold text-slate-700")
                self.job_detail = ui.label("No operation is running").classes("text-[11px] text-slate-500 truncate")
            ui.space()
            self.job_count = ui.label("").classes("text-[11px] text-slate-500")
            self.job_progress = ui.linear_progress(value=0, show_value=False, color="primary").classes("w-44")
            self.job_progress.set_visibility(False)
            with ui.row().classes("items-center gap-2 shrink-0 text-primary"):
                ui.icon("table_view", size="18px")
                ui.label("Workspace data").classes("text-xs font-semibold")
                ui.icon("expand_less", size="20px")

        with ui.column().classes("page-shell gap-0"):
            with ui.tab_panels(self.nav_tabs, value=self.tabs["load"], animated=False, keep_alive=True).classes("page-panels w-full"):
                with ui.tab_panel(self.tabs["load"]): self._build_load_page()
                with ui.tab_panel(self.tabs["preprocess"]): self._build_preprocess_page()
                with ui.tab_panel(self.tabs["web"]): self._build_web_page()
                with ui.tab_panel(self.tabs["units"]): self._build_unit_page()
                with ui.tab_panel(self.tabs["columns"]): self._build_column_page()
                with ui.tab_panel(self.tabs["deduplicate"]): self._build_deduplication_page()
                with ui.tab_panel(self.tabs["search"]): self._build_search_page()
                with ui.tab_panel(self.tabs["export"]): self._build_export_page()
        with ui.dialog().props(
            'position=bottom transition-show=slide-up transition-hide=slide-down transition-duration=280'
        ).classes("data-dialog") as self.data_dialog:
            self._build_overview()
        ui.timer(0.2, self._refresh_job_status)

    @staticmethod
    def _page_heading(title: str) -> None:
        with ui.column().classes("gap-1 mb-5"):
            ui.label(title).classes("page-title")

    def _build_load_page(self) -> None:
        self._page_heading(
            "Load a molecular dataset",
        )
        self.load_source = ui.toggle(
            {"file": "Upload file", "text": "Paste SMILES"}, value="file",
        ).props("no-caps unelevated toggle-color=primary").classes("mb-3")
        ui.label("Choose one import method. Only the selected source will be loaded.").classes("section-copy")
        with ui.column().classes("w-full"):
            with _work_card("Upload file", "CSV, Excel, TXT, SDF, MOL and SMI files are supported.") as file_card:
                file_card.bind_visibility_from(self.load_source, "value", value="file")
                self.uploader = ui.upload(
                    label="Drop a dataset here or browse", auto_upload=True, max_files=1,
                    on_upload=self._handle_upload,
                    on_rejected=lambda: ui.notify("The selected file could not be uploaded.", type="negative"),
                ).props('accept=".csv,.xlsx,.xls,.txt,.sdf,.mol,.smi" flat bordered').classes("upload-zone")
                self.upload_status = ui.label("No file selected").classes("text-xs text-slate-500")
                with ui.row().classes("items-center gap-5"):
                    self.header_check = ui.checkbox("First row contains column names", value=True, on_change=self._redetect_upload)
                    self.header_check.set_visibility(False)
                    self.sheet_select = ui.select([], label="Excel sheet", on_change=self._sheet_changed).props("outlined dense").classes("min-w-48")
                    self.sheet_select.set_visibility(False)
                ui.separator().classes("my-1")
                ui.label("Column mapping").classes("section-title")
                with ui.row().classes("grid-3"):
                    for label, key in self.MAPPINGS:
                        self.mapping[key] = ui.select([NONE_TEXT], label=label, value=NONE_TEXT).props("outlined dense options-dense").classes("w-full")
            with _work_card("Paste SMILES", "Use one structure per line for quick inspection or small ad hoc datasets.") as text_card:
                text_card.bind_visibility_from(self.load_source, "value", value="text")
                ui.label("Optional values must match the SMILES lines. Leave a value line blank for a missing value.").classes("section-copy")
                with ui.row().classes("grid-2"):
                    with ui.column().classes("w-full"):
                        self.smiles_text = ui.textarea(
                            label="SMILES list", placeholder="CC\nCCC\nCCCC",
                        ).props("outlined rows=8").classes("w-full")
                        self.text_column = ui.input("SMILES column name", value="Smiles").props("outlined dense").classes("w-full")
                    with ui.column().classes("w-full"):
                        self.values_text = ui.textarea(
                            label="Values (optional)", placeholder="1.2\n3.4\n5.6",
                        ).props("outlined rows=8").classes("w-full")
                        self.value_column = ui.input("Value column name", value="Value").props("outlined dense").classes("w-full")
            with ui.row().classes("action-row"):
                ui.button("Load dataset", icon="database", on_click=self._load_selected_source).props("unelevated no-caps color=primary").classes("primary-action")

    async def _load_selected_source(self) -> None:
        if self.load_source.value == "file":
            await self._load_uploaded_file()
        else:
            await self._load_smiles_text()

    def _build_preprocess_page(self) -> None:
        self._page_heading(
            "Standardize molecular structures",
        )
        with ui.row().classes("grid-2 items-start"):
            with _work_card("Fragment handling", "Choose how salts, solvents and disconnected components are treated."):
                with ui.row().classes("grid-2"):
                    self.remove_salts = ui.checkbox("Remove salts", value=True)
                    self.remove_solvents = ui.checkbox("Remove solvents", value=True)
                    self.remove_inorganic = ui.checkbox("Remove inorganic structures", value=True)
                self.mixture = ui.select(
                    {"reject": "Reject unresolved components", "keep": "Preserve remaining components", "largest": "Extract largest component"},
                    label="Mixture handling", value="reject",
                ).props("outlined dense").classes("w-full")
                self.hac = ui.number("Minimum heavy atoms for largest component", value=3, min=0, max=999, precision=0).props("outlined dense").classes("w-full")
                self.mixture.tooltip("All-identical components collapse to one before every mode, including Keep. Largest requires a unique largest identity; distinct ties are rejected as Ambiguous parent.")
            with _work_card("Normalization and validation", "Set charge, atom and representation policies for modeling-ready structures."):
                with ui.row().classes("grid-2"):
                    self.sanitize = ui.checkbox("RDKit sanitize", value=True)
                    self.remove_hs = ui.checkbox("Remove explicit hydrogens", value=True)
                    self.add_hs = ui.checkbox("Add explicit hydrogens at the end", value=False)
                    self.add_hs.tooltip("When both hydrogen options are enabled, remove explicit hydrogens first and add them after all chemical processing.")
                    self.remove_isotopes = ui.checkbox("Remove isotopes", value=True)
                    self.remove_stereo = ui.checkbox("Remove stereochemistry", value=False)
                    self.reject_radicals = ui.checkbox("Reject radicals", value=True)
                    self.neutralize = ui.checkbox("Neutralize charges", value=True)
                    self.reject_non_neutral = ui.checkbox("Reject non-neutral structures", value=False)
                self.element_policy = ui.select(
                    {"allow_all": "Allow all parsed elements", "reject_metals": "Reject metals only", "allowed_atoms": "Use allowed-atom list"},
                    label="Element policy", value="allow_all",
                ).props("outlined dense").classes("w-full")
        with _work_card("Execution"):
            with ui.row().classes("items-end gap-3 flex-wrap w-full"):
                self.worker_jobs = ui.number(
                    "Worker processes", value=1, min=1, max=max(1, multiprocessing.cpu_count()), precision=0,
                ).props("outlined dense").classes("w-44")
                ui.space()
                ui.button("Manage chemical rules", icon="tune", on_click=self._open_rules_dialog).props("outline no-caps color=primary").classes("primary-action")
                ui.button("Run preprocessing", icon="play_arrow", on_click=self._run_preprocessing).props("unelevated no-caps color=primary").classes("primary-action")
        self._build_rules_dialog()

    def _build_rules_dialog(self) -> None:
        with ui.dialog() as self.rules_dialog, ui.card().classes("rules-dialog-card w-[960px] max-w-[94vw] p-0"):
            with ui.row().classes("items-start w-full px-6 pt-5"):
                with ui.column().classes("gap-0"):
                    ui.label("Chemical processing rules").classes("text-xl font-bold text-slate-800")
                    ui.label("Edit the rule set used by preprocessing.").classes("text-xs text-slate-500")
                ui.space()
                ui.button(icon="close", on_click=self.rules_dialog.close).props("flat round")
            with ui.tabs().props("no-caps align=justify").classes("rules-tabs segmented-tabs") as rules_tabs:
                atom_tab = ui.tab("atoms", "Atoms")
                salt_tab = ui.tab("salts", "Salts")
                solvent_tab = ui.tab("solvents", "Solvents")
                neutral_tab = ui.tab("neutralization", "Neutralization")
            with ui.tab_panels(rules_tabs, value=atom_tab, animated=False, keep_alive=True).classes("rules-panels w-full px-6 pb-6 pt-2"):
                self.rule_inputs: dict[str, Any] = {}
                for tab, kind, label, placeholder in (
                    (atom_tab, "atoms", "Atoms", "Si, Zr"),
                    (salt_tab, "salts", "SMARTS patterns", "[Hg+2]"),
                    (solvent_tab, "solvents", "SMILES strings", "CCC"),
                ):
                    with ui.tab_panel(tab), ui.row().classes("rules-layout"):
                        with ui.column().classes("w-full gap-2"):
                            ui.label(f"Current {label.lower()}").classes("text-xs font-semibold text-slate-600")
                            current = ui.element("div").classes("rule-list")
                            if kind == "atoms":
                                current.classes("rule-list--atoms")
                        with ui.column().classes("w-full gap-3"):
                            entry = ui.input(label, placeholder=placeholder).props("outlined dense").classes("w-full")
                            self.rule_inputs[kind] = (current, entry)
                            with ui.row().classes("gap-2"):
                                ui.button("Add", on_click=lambda _=None, k=kind: self._change_simple_rule(k, True)).props("unelevated no-caps color=primary")
                                ui.button("Remove", on_click=lambda _=None, k=kind: self._change_simple_rule(k, False)).props("outline no-caps color=primary")
                with ui.tab_panel(neutral_tab), ui.row().classes("rules-layout"):
                    with ui.column().classes("w-full gap-2"):
                        ui.label("Current neutralization rules").classes("text-xs font-semibold text-slate-600")
                        self.neutral_rules = ui.element("div").classes("rule-list rule-list--neutral")
                    with ui.column().classes("w-full gap-3"):
                        self.neutral_reactant = ui.input("Reactant SMARTS").props("outlined dense").classes("w-full")
                        self.neutral_product = ui.input("Product SMILES").props("outlined dense").classes("w-full")
                        with ui.row().classes("gap-2"):
                            ui.button("Add", on_click=lambda: self._change_neutral_rule(True)).props("unelevated no-caps color=primary")
                            ui.button("Remove by reactant", on_click=lambda: self._change_neutral_rule(False)).props("outline no-caps color=primary")

    def _build_web_page(self) -> None:
        self._page_heading(
            "Request identifiers and properties",
        )
        with _work_card("Request plan", "Select at least one source, requested property and input identifier."):
            with ui.row().classes("grid-3"):
                self.web_sources = ui.select(
                    ["pubchem", "chemspider", "comptox", "cas", "cactus", "chembl"], label="Data sources",
                    value=["pubchem"], multiple=True, with_input=True, clearable=True,
                ).props("outlined dense use-chips options-dense").classes("w-full")
                self.web_outputs = ui.select(
                    ["smiles", "cas", "name", "iupac", "mw"], label="Requested properties",
                    value=["cas", "iupac"], multiple=True, with_input=True, clearable=True,
                ).props("outlined dense use-chips options-dense").classes("w-full")
                self.web_inputs = ui.select(
                    ["smiles", "cas", "name"], label="Input identifiers", value=["smiles"], multiple=True, with_input=True, clearable=True,
                ).props("outlined dense use-chips options-dense").classes("w-full")
        with _work_card("Connection settings", "Defaults are conservative. Increase concurrency only when the selected service permits it."):
            with ui.row().classes("settings-grid"):
                self.web_interval = ui.number("Request interval", value=0.3, min=0, max=3600, step=0.1, suffix=" s").props("outlined dense")
                self.web_retries = ui.number("Retries", value=3, min=0, max=100, precision=0).props("outlined dense")
                self.web_delay = ui.number("Retry delay", value=30, min=0, max=86400, precision=0, suffix=" s").props("outlined dense")
                self.web_workers = ui.number("Concurrent requests", value=4, min=1, max=64, precision=0).props("outlined dense")
                self.web_batch = ui.number("Batch limit", value=1500, min=0, max=10_000_000, precision=0).props("outlined dense")
                self.web_rest = ui.number("Rest duration", value=300, min=0, max=86400, precision=0, suffix=" s").props("outlined dense")
            with ui.expansion("API keys", caption="Optional credentials for ChemSpider, CompTox and CAS", icon="key").classes("w-full"):
                with ui.row().classes("grid-3 p-2"):
                    self.chemspider_key = ui.input("ChemSpider API key", password=True, password_toggle_button=True).props("outlined dense")
                    self.comptox_key = ui.input("CompTox API key", password=True, password_toggle_button=True).props("outlined dense")
                    self.cas_key = ui.input("CAS API key", password=True, password_toggle_button=True).props("outlined dense")
                self.force_api = ui.checkbox("Force raw API mode", value=False).classes("px-2")
            with ui.row().classes("action-row"):
                ui.button("Start web request", icon="cloud_download", on_click=self._run_web_request).props("unelevated no-caps color=primary").classes("primary-action")

    def _build_unit_page(self) -> None:
        self._page_heading(
            "Standardize measurement units",
        )
        self.unit_context = ui.label("Load data with target and unit columns to configure this step.").classes("context-strip w-full mb-4")
        with _work_card('Columns to transform', 'For condition columns, choose a value/unit pair. Conversion runs before base-10 log transformation; original columns are preserved.'):
            with ui.row().classes('grid-2'):
                self.unit_scope = ui.select(['Primary target', 'Condition column'], value='Primary target', label='Apply to', on_change=lambda _: self._refresh_data_controls()).props('outlined dense').classes('w-full')
                self.condition_log = ui.select(['None', 'log10', '-log10'], value='None', label='Condition transformation').props('outlined dense').classes('w-full')
                self.condition_value = ui.select([], label='Value column').props('outlined dense').classes('w-full')
                self.condition_unit = ui.select([], label='Unit column', on_change=lambda _: self.unit_rules_refresh.refresh()).props('outlined dense').classes('w-full')
            ui.label('Leave Standard unit empty to apply only log10 / -log10. Invalid rows enter the exclusion table.').classes('section-copy')
        with _work_card("Conversion target", "Formulas use x for the source value and mw for molecular weight."):
            with ui.row().classes("grid-3"):
                self.unit_standard = ui.input(
                    label="Standard unit", placeholder="Enter a unit, e.g. mol/L", autocomplete=[],
                    on_change=self._unit_target_changed,
                ).props("outlined dense").classes("w-full")
                self.mw_basis = ui.select(
                    {"original": "Original input structure", "standardized": "Standardized parent", "column": "Dataset MW column"},
                    label="Molecular-weight basis", value=None, on_change=self._mw_basis_changed,
                ).props("outlined dense").classes("w-full")
                self.mw_basis.tooltip("Choose a basis for every mass/molar conversion. Standardized parent requires preprocessing; conversions without molecular weight can leave this blank.")
                self.mw_column = ui.select([], label="Molecular-weight column").props("outlined dense").classes("w-full")
                self.mw_column.set_enabled(False)
            ui.separator().classes("my-2")

            @ui.refreshable
            def unit_rules() -> None:
                target = _optional_text(self.unit_standard.value)
                self.unit_rule_inputs = {}
                sources = [unit for unit in self._detected_units() if target and unit != target]
                if not sources:
                    ui.label("No conversion formulas are required for the current target.").classes("text-sm text-slate-500 py-4")
                    return
                with ui.column().classes("w-full gap-2"):
                    with ui.row().classes("w-full text-[11px] font-semibold text-slate-500 px-1"):
                        ui.label("Source").classes("w-1/5")
                        ui.label("Target").classes("w-1/5")
                        ui.label("Formula").classes("flex-1")
                    processor = UnitProcessor()
                    for source in sources:
                        key = (source, target)
                        formula = self.state.unit_rules.get(key, processor.get_rule(source, target) or "")
                        with ui.row().classes("w-full items-center gap-3"):
                            ui.label(source).classes("w-1/5 text-sm font-medium")
                            ui.label(target).classes("w-1/5 text-sm text-slate-600")
                            entry = ui.input(value=formula, placeholder="x * 1000").props("outlined dense").classes("flex-1")
                            self.unit_rule_inputs[key] = entry

            self.unit_rules_refresh = unit_rules
            unit_rules()
            with ui.row().classes("action-row"):
                self.unit_run_button = ui.button(
                    "Run unit standardization", icon="swap_horiz", on_click=self._run_unit_standardization,
                ).props("unelevated no-caps color=primary").classes("primary-action")

    def _build_column_page(self) -> None:
        self._page_heading('Adjust column values before deduplication')
        self._build_value_filter()
        self._build_value_merge()

    def _build_deduplication_page(self) -> None:
        self._page_heading(
            "Create one modeling record per structure",
        )
        with _work_card("Grouping", "Condition columns define distinct experimental contexts and are not merged across groups."):
            self.dedup_conditions = ui.select(
                [], label="Condition columns", value=[], multiple=True, clearable=True, with_input=True,
            ).props("outlined dense use-chips options-dense").classes("w-full")
            self.dedup_drop_na = ui.checkbox("Drop records with missing condition values", value=False)
        with _work_card("Resolution method", "Discrete targets support majority voting or explicit priority selection."):
            with ui.row().classes("grid-2"):
                self.dedup_type = ui.select(
                    ["continuous", "discrete", "smiles"], label="Data type", value="continuous", on_change=self._dedup_type_changed,
                ).props("outlined dense").classes("w-full")
                self.dedup_method = ui.select(["auto", "3sigma", "IQR"], label="Outlier filtering", value="auto").props("outlined dense").classes("w-full")
                self.dedup_method.tooltip("Continuous targets: transform first, filter outliers, then aggregate. Groups of 3 or fewer skip filtering.")
                self.dedup_aggregation = ui.select(["mean", "max", "min"], label="Aggregation", value="mean").props("outlined dense").classes("w-full")
                self.dedup_transform = ui.select(["None", "-log10", "log10"], label="Target transformation", value="None").props("outlined dense").classes("w-full")
                self.dedup_priority = ui.input("Priority values", placeholder="Active, Intermediate").props("outlined dense").classes("w-full")
                self.dedup_priority.bind_enabled_from(self.dedup_method, "value", backward=lambda value: value == "priority")
                self.dedup_priority.tooltip("Used only by priority, in the listed order. If none match, fall back to voting.")
            with ui.row().classes("action-row"):
                ui.button("Run deduplication", icon="call_merge", on_click=self._run_deduplication).props("unelevated no-caps color=primary").classes("primary-action")

    def _build_search_page(self) -> None:
        self._page_heading(
            "Search and filter records",
        )
        with ui.row().classes("grid-2 items-start"):
            with _work_card("Substructure search", "The result is added as a Boolean column and the source records are retained."):
                self.search_pattern = ui.input("SMARTS or SMILES pattern", placeholder="c1ccccc1").props("outlined").classes("w-full")
                self.search_smarts = ui.checkbox("Interpret pattern as SMARTS", value=True)
                with ui.row().classes("action-row"):
                    ui.button("Run search", icon="search", on_click=self._run_substructure_search).props("unelevated no-caps color=primary").classes("primary-action")
            with _work_card("Atom-count filter", "Enabled ranges are inclusive. Records outside a range are removed from the working dataset."):
                self.use_heavy = ui.checkbox("Apply heavy-atom range", value=True)
                with ui.row().classes("grid-2"):
                    self.min_heavy = ui.number("Minimum heavy atoms", value=0, min=0, max=999999, precision=0).props("outlined dense")
                    self.max_heavy = ui.number("Maximum heavy atoms", value=999, min=0, max=999999, precision=0).props("outlined dense")
                self.use_total = ui.checkbox("Apply total-atom range", value=False)
                with ui.row().classes("grid-2"):
                    self.min_total = ui.number("Minimum total atoms", value=0, min=0, max=999999, precision=0).props("outlined dense")
                    self.max_total = ui.number("Maximum total atoms", value=999, min=0, max=999999, precision=0).props("outlined dense")
                with ui.row().classes("action-row"):
                    ui.button("Apply filter", icon="filter_alt", on_click=self._run_atom_filter).props("unelevated no-caps color=primary").classes("primary-action")

    def _build_value_filter(self) -> None:
        with _work_card("Filter by column values", "Select a column to see every distinct value and its row count. No selected values keeps all rows; removed rows enter the exclusion table."):
            self.value_filter_column = ui.select([], label="Column", on_change=self._value_filter_column_changed).props("outlined dense").classes("w-full")
            self.value_filter_mode = ui.radio({'keep': 'Keep selected values', 'remove': 'Remove selected values'}, value='keep').props('inline')
            self.value_filter_values = ui.select([], label="Values", value=[], multiple=True, with_input=True, clearable=True).props("outlined dense use-chips").classes("w-full")
            self.value_filter_summary = ui.label("Choose a column.").classes("section-copy")
            ui.button("Apply column filter", icon="filter_alt", on_click=self._run_value_filter).props("unelevated no-caps color=primary")

    def _build_value_merge(self) -> None:
        with _work_card('Merge column values', 'Combine selected values into one label for grouping. Original columns and all rows are preserved; select the new column under Condition columns.'):
            self.value_merge_column = ui.select([], label='Column', on_change=self._value_merge_column_changed).props('outlined dense').classes('w-full')
            self.value_merge_mode = ui.radio({'single': 'One group', 'batch': 'Multiple groups (JSON)'}, value='single').props('inline')
            self.value_merge_values = ui.select([], label='Values to merge', value=[], multiple=True, with_input=True, clearable=True).props('outlined dense use-chips').classes('w-full')
            self.value_merge_values.bind_visibility_from(self.value_merge_mode, 'value', value='single')
            self.value_merge_groups = ui.textarea('Merge rules (JSON array)', placeholder='[{"values": ["Embryo", "Egg"], "replacement": "Embryonic"},\n {"values": ["Larva", "Fry"], "replacement": "Larval / post hatch"}]').props('outlined autogrow input-style="font-family:monospace"').classes('w-full')
            self.value_merge_groups.bind_visibility_from(self.value_merge_mode, 'value', value='batch')
            with ui.row().classes('grid-2'):
                self.value_merge_label = ui.input('Replace with', value='other').props('outlined dense').classes('w-full')
                self.value_merge_label.bind_visibility_from(self.value_merge_mode, 'value', value='single')
                self.value_merge_output = ui.input('New column name (optional)', placeholder='Source column (Merged)').props('outlined dense').classes('w-full')
            self.value_merge_summary = ui.label('Choose a column. Empty selection makes no changes.').classes('section-copy')
            ui.label('Batch rules match the original values simultaneously. Unlisted values stay unchanged; conflicting groups are rejected.').classes('section-copy').bind_visibility_from(self.value_merge_mode, 'value', value='batch')
            ui.button('Apply merge', icon='merge_type', on_click=self._run_value_merge).props('unelevated no-caps color=primary')

    def _value_merge_column_changed(self, _event=None) -> None:
        column = self.value_merge_column.value
        self._merge_value_options = (self.pipeline.get_column_values(column)
                                     if self.pipeline.df is not None and column in self.pipeline.df else [])
        options = {index: f"{'(missing)' if entry['value'] is None else repr(entry['value'])}  ({entry['count']:,} rows)"
                   for index, entry in enumerate(self._merge_value_options)}
        self.value_merge_values.set_options(options, value=[])
        self.value_merge_summary.set_text(f'{len(options):,} distinct values. Empty selection makes no changes.')

    async def _run_value_merge(self) -> None:
        if not self._require_data():
            return
        column = self.value_merge_column.value
        if column not in self.pipeline.df:
            ui.notify('Choose a column.', type='warning')
            return
        values = [self._merge_value_options[index]['value'] for index in (self.value_merge_values.value or [])]
        groups = None
        if self.value_merge_mode.value == 'batch':
            try:
                groups = json.loads(self.value_merge_groups.value or '[]')
                if not isinstance(groups, list):
                    raise ValueError('Rules must be a JSON array.')
                from .value_filter import merge_rules
                merge_rules(groups=groups)
            except (ValueError, TypeError) as exc:
                ui.notify(f'Invalid merge rules: {exc}', type='warning')
                return
            values = []
        if not values and not groups:
            ui.notify('No values selected; no changes made.', type='info')
            return
        replacement = self.value_merge_label.value if groups is None else 'other'
        if not replacement or not replacement.strip():
            ui.notify('Enter a replacement label, for example other.', type='warning')
            return
        output = _optional_text(self.value_merge_output.value) or column + ' (Merged)'

        def action(pipeline, _progress, _status):
            pipeline.merge_column_values(column, values, replacement, output, groups=groups)

        await self._run_pipeline_job('Merge column values', action,
                                     lambda pipeline: f'Created {output}; select it as a deduplication condition.')

    def _value_filter_column_changed(self, _event=None) -> None:
        column = self.value_filter_column.value
        if self.pipeline.df is None or column not in self.pipeline.df.columns:
            self._column_value_options = []
        else:
            self._column_value_options = self.pipeline.get_column_values(column)
        options = {
            index: f"{'(missing)' if entry['value'] is None else repr(entry['value'])}  ({entry['count']:,} rows)"
            for index, entry in enumerate(self._column_value_options)
        }
        self.value_filter_values.set_options(options, value=[])
        self.value_filter_summary.set_text(f"{len(options):,} distinct values. Empty selection keeps all rows.")

    async def _run_value_filter(self) -> None:
        if not self._require_data():
            return
        column = self.value_filter_column.value
        if column not in self.pipeline.df.columns:
            ui.notify("Choose a column.", type="warning")
            return
        values = [self._column_value_options[index]['value'] for index in (self.value_filter_values.value or [])]
        if not values:
            ui.notify("No values selected; all rows are retained.", type="info")
            return
        mode = self.value_filter_mode.value

        def action(pipeline, _progress, _status):
            pipeline.filter_by_values(column, values, mode)

        await self._run_pipeline_job("Column value filter", action,
                                     lambda pipeline: f"Filter complete: {len(pipeline.df):,} records remain.")

    def _build_export_page(self) -> None:
        self._page_heading(
            "Export modeling-ready data",
        )
        with _work_card("Processed dataset", "Recommended columns preserve identifiers, standardized structures, targets and audit fields."):
            with ui.row().classes("items-end gap-4 w-full"):
                self.export_filename = ui.input("File name", value="diptox-processed", placeholder="My dataset").props("outlined dense").classes("min-w-64")
                self.export_filename.tooltip("The selected format determines the extension. Excluded records use the same name with -excluded.csv.")
                self.export_format = ui.select(["csv", "xlsx", "txt", "sdf", "smi"], label="Format", value="csv").props("outlined dense").classes("w-40")
                with ui.row().classes("gap-1"):
                    ui.button("Recommended", on_click=self._select_recommended_columns).props("flat no-caps color=primary")
                    ui.button("Select all", on_click=self._select_all_columns).props("flat no-caps color=primary")
                    ui.button("Clear", on_click=self._clear_export_columns).props("flat no-caps color=secondary")
            self.export_columns = ui.select([], label="Columns", value=[], multiple=True, with_input=True, clearable=True, on_change=self._export_columns_changed).props("outlined dense use-chips options-dense").classes("w-full")
            with ui.row().classes("action-row"):
                ui.button("Download results", icon="download", on_click=self._export_results).props("unelevated no-caps color=primary").classes("primary-action")
        with _work_card("Excluded records", "Includes invalid structures from preprocessing and records excluded during deduplication."):
            self.excluded_summary = ui.label("No excluded records are available.").classes("text-sm text-slate-500")
            with ui.row().classes("action-row"):
                self.export_excluded_button = ui.button(
                    "Download excluded records", icon="file_download", on_click=self._export_excluded,
                ).props("outline no-caps color=primary")

    def _build_overview(self) -> None:
        @ui.refreshable
        def overview() -> None:
            frame = self.pipeline.df
            with ui.card().classes("work-card data-surface w-full"):
                with ui.row().classes("items-center w-full"):
                    with ui.column().classes("gap-0"):
                        ui.label("Workspace data").classes("section-title")
                        summary = "No dataset loaded" if frame is None else f"{len(frame):,} rows · {len(frame.columns):,} columns · preview limited to 200 rows"
                        ui.label(summary).classes("section-copy")
                    ui.space()
                    if frame is not None:
                        key = "Canonical SMILES" if self.pipeline._preprocess_key else (self.pipeline.smiles_col or "Not mapped")
                        ui.label(f"Structure field: {key}").classes("text-xs text-slate-500")
                    ui.button(icon="close", on_click=self.data_dialog.close).props('flat round aria-label="Close workspace data"')
                if frame is None:
                    with ui.column().classes("empty-state w-full"):
                        ui.icon("table_view", size="34px").classes("text-slate-400")
                        ui.label("Load a dataset to inspect records and processing history.").classes("text-sm mt-2")
                    return
                with ui.tabs().props("no-caps align=justify").classes("segmented-tabs w-full") as data_tabs:
                    preview_tab = ui.tab("preview", "Dataset")
                    history_tab = ui.tab("history", "Processing history")
                    excluded_tab = ui.tab("excluded", f"Excluded ({len(self.pipeline.excluded_df):,})")
                with ui.tab_panels(data_tabs, value=preview_tab, animated=False).classes("w-full"):
                    with ui.tab_panel(preview_tab):
                        ui.aggrid(_grid_options(frame), theme="quartz").classes("w-full").style("height:360px")
                    with ui.tab_panel(history_tab):
                        history = self.pipeline.get_history()
                        if history.empty:
                            ui.label("No processing steps have been recorded yet.").classes("text-sm text-slate-500 py-8")
                        else:
                            ui.aggrid(_grid_options(history, 100), theme="quartz").classes("w-full").style("height:300px")
                    with ui.tab_panel(excluded_tab):
                        if self.pipeline.excluded_df.empty:
                            ui.label("No invalid or excluded records are available.").classes("text-sm text-slate-500 py-8")
                        else:
                            ui.aggrid(_grid_options(self.pipeline.excluded_df), theme="quartz").classes("w-full").style("height:300px")
        self.overview_refresh = overview
        overview()

    def _navigation_changed(self, event: events.ValueChangeEventArguments) -> None:
        value = event.value.name if hasattr(event.value, "name") else str(event.value)
        label = next((item[1] for item in PAGE_DEFINITIONS if item[0] == value), "Workflow")
        self.header_stage.set_text(label)

    def _refresh_job_status(self) -> None:
        busy, description, status, current, total, _error = self.store.snapshot()
        self.job_title.set_text(description if busy else "Ready")
        self.job_detail.set_text(status)
        if busy:
            self.status_dot.set_content('<div class="status-dot busy"></div>')
            self.job_progress.set_visibility(True)
            if total > 0:
                self.job_progress.set_value(min(max(current / total, 0.0), 1.0))
                self.job_count.set_text(f"{current:,} / {total:,}")
            else:
                self.job_progress.set_value(0.08)
                self.job_count.set_text("Working")
        else:
            self.status_dot.set_content('<div class="status-dot"></div>')
            self.job_progress.set_visibility(False)
            self.job_count.set_text("")
        self.undo_button.set_enabled(not busy and bool(self.pipeline._history))

    def _require_data(self) -> bool:
        if self.pipeline.df is None:
            ui.notify("Load a dataset first.", type="warning")
            return False
        return True

    async def _run_pipeline_job(
        self, description: str, action: Callable, success: str | Callable[[DiptoxPipeline], str],
    ) -> bool:
        try:
            result = await self.store.run_pipeline_job(description, action)
        except BusyError as exc:
            ui.notify(str(exc), type="warning")
            return False
        except Exception as exc:
            ui.notify(str(exc) or "The operation failed.", type="negative", multi_line=True, timeout=0, close_button=True)
            return False
        self._refresh_data_controls()
        message = success(result) if callable(success) else success
        ui.notify(message, type="positive")
        return True

    async def _run_read_job(self, description: str, operation: Callable) -> Any:
        try:
            return await self.store.run_read_job(description, operation)
        except BusyError as exc:
            ui.notify(str(exc), type="warning")
        except Exception as exc:
            ui.notify(str(exc) or "The operation failed.", type="negative", multi_line=True, timeout=0, close_button=True)
        return None

    async def _handle_upload(self, event: events.UploadEventArguments) -> None:
        # NiceGUI 3 uses an async FileUpload; NiceGUI 2 supplies a binary stream.
        uploaded_file = getattr(event, "file", None)
        name = uploaded_file.name if uploaded_file is not None else event.name
        suffix = Path(name).suffix.lower()
        if suffix not in SUPPORTED_SUFFIXES:
            ui.notify(f"Unsupported file type: {suffix or 'unknown'}", type="negative")
            return
        upload_dir = Path(tempfile.gettempdir()) / "diptox-uploads"
        upload_dir.mkdir(parents=True, exist_ok=True)
        destination = upload_dir / f"{uuid.uuid4().hex}{suffix}"
        if uploaded_file is not None:
            await uploaded_file.save(destination)
        else:
            def save_legacy_upload() -> None:
                event.content.seek(0)
                with destination.open("wb") as output:
                    shutil.copyfileobj(event.content, output)

            await _to_thread(save_legacy_upload)
        self.state.upload_path = destination
        self.state.upload_name = name
        self.state.sheet_name = None
        self.upload_status.set_text(f"Selected: {name} · {destination.stat().st_size / 1024:.1f} KB")
        await self._detect_columns()

    async def _redetect_upload(self, _event: Any = None) -> None:
        if self.state.upload_path:
            await self._detect_columns(keep_sheets=True)

    async def _sheet_changed(self, event: events.ValueChangeEventArguments) -> None:
        if event.value and str(event.value) != self.state.sheet_name:
            self.state.sheet_name = str(event.value)
            await self._detect_columns(keep_sheets=True)

    @staticmethod
    def _inspect_file(
        path: Path, header: bool, sheet_name: Optional[str],
    ) -> tuple[list[str], list[str], Optional[str]]:
        suffix = path.suffix.lower()
        sheets: list[str] = []
        selected_sheet = sheet_name
        if suffix in {".xlsx", ".xls"}:
            with pd.ExcelFile(path) as excel:
                sheets = list(excel.sheet_names)
                if selected_sheet not in sheets:
                    selected_sheet = sheets[0] if sheets else None
                frame = pd.read_excel(path, sheet_name=selected_sheet, nrows=5, header=0)
                columns = [str(column) for column in frame.columns]
        elif suffix in {".sdf", ".mol"}:
            from rdkit import Chem
            if suffix == ".mol":
                molecule = Chem.MolFromMolFile(str(path))
            else:
                with path.open("rb") as handle:
                    molecule = next(Chem.ForwardSDMolSupplier(handle), None)
            columns = list(map(str, molecule.GetPropsAsDict())) if molecule else []
        else:
            separator = "\t" if suffix in {".txt", ".smi"} else ","
            frame = pd.read_csv(path, sep=separator, nrows=5, header=0 if header else None)
            columns = [str(column) for column in frame.columns]
        return columns, sheets, selected_sheet

    async def _detect_columns(self, keep_sheets: bool = False) -> None:
        if not self.state.upload_path:
            return
        self.header_check.set_visibility(self.state.upload_path.suffix.lower() in HEADER_OPTION_SUFFIXES)
        try:
            columns, sheets, selected_sheet = await _to_thread(
                self._inspect_file, self.state.upload_path, bool(self.header_check.value), self.state.sheet_name,
            )
        except Exception as exc:
            self.state.detected_columns = []
            self._set_mapping_options([])
            ui.notify(f"Could not read the file header: {exc}", type="negative")
            return
        self.state.detected_columns = columns
        self.state.sheet_name = selected_sheet
        if sheets:
            self.sheet_select.set_options(sheets, value=selected_sheet)
            self.sheet_select.set_visibility(len(sheets) > 1)
        elif not keep_sheets:
            self.sheet_select.set_visibility(False)
        self._set_mapping_options(columns)
        self.upload_status.set_text(f"Selected: {self.state.upload_name} · {len(columns)} column(s) detected")

    def _set_mapping_options(self, columns: Iterable[str]) -> None:
        columns = list(columns)
        hints = {
            "smiles": ("smiles", "canonical smiles", "canonical_smiles"),
            "name": ("name", "chemical name", "compound name"),
            "cas": ("cas", "casrn", "cas number"),
            "id": ("id", "compound id"),
            "target": ("target", "value", "activity"),
            "unit": ("unit", "units"),
        }
        lower = {column.lower(): column for column in columns}
        options = [NONE_TEXT, *columns]
        for key, select in self.mapping.items():
            old = select.value
            preferred = old if old in columns else next((lower[name] for name in hints[key] if name in lower), NONE_TEXT)
            select.set_options(options, value=preferred)

    async def _load_uploaded_file(self) -> None:
        path = self.state.upload_path
        if path is None:
            ui.notify("Upload a dataset first.", type="warning")
            return
        values = {key: _optional_text(select.value) for key, select in self.mapping.items()}
        if path.suffix.lower() in {".sdf", ".mol"} and values["smiles"] is None:
            values["smiles"] = "smiles"
        chosen = [value for value in values.values() if value]
        if len(chosen) != len(set(chosen)):
            ui.notify("A source column is assigned to more than one role.", type="negative")
            return
        kwargs: dict[str, Any] = {}
        suffix = path.suffix.lower()
        if suffix in HEADER_OPTION_SUFFIXES:
            kwargs["header"] = 0 if self.header_check.value else None
        if suffix in {".csv", ".txt"} and not self.header_check.value:
            kwargs["names"] = self.state.detected_columns
        if suffix in {".xlsx", ".xls"} and self.state.sheet_name:
            kwargs["sheet_name"] = self.state.sheet_name

        def action(pipeline: DiptoxPipeline, _progress: Callable, status: Callable) -> None:
            status("Reading dataset")
            pipeline.load_data(
                input_data=str(path), smiles_col=values["smiles"], name_col=values["name"],
                cas_col=values["cas"], id_col=values["id"], target_col=values["target"],
                unit_col=values["unit"], **kwargs,
            )
        await self._run_pipeline_job("Loading data", action, lambda pipeline: f"Loaded {len(pipeline.df):,} records.")

    async def _load_smiles_text(self) -> None:
        try:
            frame, column, target = self._parse_pasted_data(
                self.smiles_text.value or "", self.text_column.value or "",
                self.values_text.value or "", self.value_column.value or "",
            )
        except ValueError as error:
            ui.notify(str(error), type="warning")
            return

        def action(pipeline: DiptoxPipeline, _progress: Callable, status: Callable) -> None:
            status("Loading SMILES list")
            pipeline.load_data(frame, smiles_col=column, target_col=target)
        await self._run_pipeline_job("Loading data", action, lambda pipeline: f"Loaded {len(pipeline.df):,} records.")

    @staticmethod
    def _parse_pasted_data(smiles: str, column: str, values: str, value_column: str):
        column, value_column = column.strip(), value_column.strip()
        lines = smiles.splitlines()
        if not any(line.strip() for line in lines) or not column:
            raise ValueError("Enter at least one SMILES string and a column name.")
        if not values.strip():
            return pd.DataFrame({column: [line.strip() for line in lines if line.strip()]}), column, None
        if not value_column or value_column == column:
            raise ValueError("Use a distinct, non-empty value column name.")
        value_lines = values.splitlines()
        if len(lines) != len(value_lines):
            raise ValueError(f"Line counts must match: {len(lines)} SMILES lines and {len(value_lines)} value lines.")
        records = []
        for index, (smiles_line, value_line) in enumerate(zip(lines, value_lines), 1):
            if not smiles_line.strip():
                if value_line.strip():
                    raise ValueError(f"Line {index} has a value but no SMILES.")
                continue
            records.append({column: smiles_line.strip(), value_column: value_line.strip() or None})
        frame = pd.DataFrame(records)
        # Preserve categorical targets; convert the column only if all supplied values are numeric.
        numeric = pd.to_numeric(frame[value_column], errors="coerce")
        if numeric.notna().sum() == frame[value_column].notna().sum():
            frame[value_column] = numeric
        return frame, column, value_column

    async def _run_preprocessing(self) -> None:
        if not self._require_data():
            return
        config = {
            "remove_salts": bool(self.remove_salts.value), "remove_solvents": bool(self.remove_solvents.value),
            "remove_inorganic": bool(self.remove_inorganic.value), "mixture_mode": self.mixture.value,
            "hac_threshold": int(self.hac.value or 0) if self.mixture.value == "largest" else 0,
            "sanitize": bool(self.sanitize.value), "remove_hs": bool(self.remove_hs.value),
            "add_hs": bool(self.add_hs.value),
            "remove_isotopes": bool(self.remove_isotopes.value), "remove_stereo": bool(self.remove_stereo.value),
            "reject_radical_species": bool(self.reject_radicals.value), "neutralize": bool(self.neutralize.value),
            "reject_non_neutral": bool(self.reject_non_neutral.value), "element_policy": self.element_policy.value,
            "n_jobs": int(self.worker_jobs.value or 1),
        }

        def action(pipeline: DiptoxPipeline, progress: Callable, status: Callable) -> None:
            status("Processing molecules")
            pipeline.preprocess(progress_callback=progress, **config)
        await self._run_pipeline_job("Preprocessing", action, lambda pipeline: f"Preprocessing complete: {len(pipeline.df):,} records.")

    def _open_rules_dialog(self) -> None:
        if self.store.busy:
            ui.notify("Wait for the current operation before editing rules.", type="warning")
            return
        self._refresh_rule_dialog()
        self.rules_dialog.open()

    def _refresh_rule_dialog(self) -> None:
        rules = self.pipeline.chem_processor.get_current_rules_dict()
        for kind in ("atoms", "salts", "solvents"):
            current, _entry = self.rule_inputs[kind]
            current.clear()
            with current:
                if not rules[kind]:
                    ui.label("No rules").classes("rule-empty")
                for rule in rules[kind]:
                    ui.label(str(rule)).classes("rule-tile")
        self.neutral_rules.clear()
        with self.neutral_rules:
            neutralization = rules["neutralization"]
            if not neutralization:
                ui.label("No neutralization rules").classes("rule-empty")
            for reactant, product in neutralization:
                ui.label(f"{reactant}  →  {product}").classes("rule-tile")

    async def _change_simple_rule(self, kind: str, add: bool) -> None:
        _current, entry = self.rule_inputs[kind]
        values = [part.strip() for part in (entry.value or "").split(",") if part.strip()]
        if not values:
            ui.notify("Enter at least one rule.", type="warning")
            return
        names = {"atoms": "manage_atom_rules", "salts": "manage_default_salt", "solvents": "manage_default_solvent"}

        def action(pipeline: DiptoxPipeline, _progress: Callable, status: Callable) -> None:
            status("Updating chemical rules")
            failed = getattr(pipeline, names[kind])(values, add=add)
            if failed:
                raise ValueError(f"Could not {'add' if add else 'remove'}: {failed}")
        success = await self._run_pipeline_job("Updating rules", action, "Chemical rules updated.")
        if success:
            entry.set_value("")
            self._refresh_rule_dialog()

    async def _change_neutral_rule(self, add: bool) -> None:
        reactant = (self.neutral_reactant.value or "").strip()
        product = (self.neutral_product.value or "").strip()
        if not reactant or (add and not product):
            ui.notify("Enter the required rule fields.", type="warning")
            return

        def action(pipeline: DiptoxPipeline, _progress: Callable, status: Callable) -> None:
            status("Updating neutralization rules")
            result = pipeline.add_neutralization_rule(reactant, product) if add else pipeline.remove_neutralization_rule(reactant)
            if not result:
                raise ValueError("The rule is invalid or was not found.")
        success = await self._run_pipeline_job("Updating rules", action, "Neutralization rules updated.")
        if success:
            self.neutral_reactant.set_value("")
            self.neutral_product.set_value("")
            self._refresh_rule_dialog()

    async def _run_web_request(self) -> None:
        if not self._require_data():
            return
        sources = list(self.web_sources.value or [])
        outputs = list(self.web_outputs.value or [])
        inputs = list(self.web_inputs.value or [])
        if not sources or not outputs or not inputs:
            ui.notify("Select at least one source, property and input identifier.", type="warning")
            return
        config = {
            "sources": sources, "chemspider_api_key": _optional_text(self.chemspider_key.value),
            "comptox_api_key": _optional_text(self.comptox_key.value), "cas_api_key": _optional_text(self.cas_key.value),
            "max_workers": int(self.web_workers.value or 1), "interval": float(self.web_interval.value or 0),
            "retries": int(self.web_retries.value or 0), "delay": int(self.web_delay.value or 0),
            "batch_limit": int(self.web_batch.value or 0), "rest_duration": int(self.web_rest.value or 0),
            "force_api_mode": bool(self.force_api.value),
        }

        def action(pipeline: DiptoxPipeline, progress: Callable, status: Callable) -> None:
            pipeline.config_web_request(status_callback=status, **config)
            pipeline.web_request(send=inputs, request=outputs, progress_callback=progress)
        await self._run_pipeline_job("Web request", action, lambda pipeline: f"Web request complete for {len(pipeline.df):,} records.")

    def _detected_units(self) -> list[str]:
        pipeline = self.pipeline
        column = self.condition_unit.value if self.unit_scope.value == 'Condition column' else pipeline.unit_col
        if pipeline.df is None or not column or column not in pipeline.df.columns:
            return []
        return [str(value) for value in pipeline.df[column].dropna().unique() if str(value).strip()]

    def _capture_unit_rules(self) -> None:
        for key, entry in self.unit_rule_inputs.items():
            formula = (entry.value or "").strip()
            if formula:
                self.state.unit_rules[key] = formula
            else:
                self.state.unit_rules.pop(key, None)

    def _unit_target_changed(self, _event: Any = None) -> None:
        self._capture_unit_rules()
        self.unit_rules_refresh.refresh()

    def _mw_basis_changed(self, _event: Any = None) -> None:
        self.mw_column.set_enabled(self.mw_basis.value == "column")

    async def _run_unit_standardization(self) -> None:
        if not self._require_data():
            return
        standard = _optional_text(self.unit_standard.value)
        condition = self.unit_scope.value == 'Condition column'
        if condition and (not self.condition_value.value or not self.condition_unit.value):
            ui.notify('Choose a value column and its unit column.', type='warning')
            return
        if not standard and not (condition and self.condition_log.value != 'None'):
            ui.notify("Choose or enter a standard unit.", type="warning")
            return
        self._capture_unit_rules()
        basis = self.mw_basis.value
        mw_column = self.mw_column.value if basis == "column" else None
        if basis == "column" and not mw_column:
            ui.notify("Choose a molecular-weight column.", type="warning")
            return
        mw_source = None if basis == "column" else basis
        rules = dict(self.state.unit_rules)

        def action(pipeline: DiptoxPipeline, _progress: Callable, status: Callable) -> None:
            status("Applying unit conversions")
            operation = pipeline.transform_column if condition else pipeline.standardize_units
            extra = dict(value_col=self.condition_value.value, unit_col=self.condition_unit.value,
                         log_transform=self.condition_log.value) if condition else {}
            operation(
                standard_unit=standard, conversion_rules=rules, molecular_weight_source=mw_source,
                molecular_weight_col=mw_column, **extra,
            )
        await self._run_pipeline_job(
            "Unit standardization", action,
            lambda pipeline: 'Condition columns created; select them under deduplication conditions.' if condition else f"Created standardized target column: {pipeline.target_col}",
        )

    def _dedup_type_changed(self, _event: Any = None) -> None:
        data_type = self.dedup_type.value
        methods = ["auto", "3sigma", "IQR"] if data_type == "continuous" else (["vote", "priority"] if data_type == "discrete" else ["auto"])
        self.dedup_method.set_options(methods, value=methods[0])
        self.dedup_method.set_label("Outlier filtering" if data_type == "continuous" else "Selection method")
        self.dedup_aggregation.set_enabled(data_type == "continuous")
        self.dedup_transform.set_enabled(data_type == "continuous")

    async def _run_deduplication(self) -> None:
        if not self._require_data():
            return
        priorities = [part.strip() for part in (self.dedup_priority.value or "").split(",") if part.strip()] or None
        if self.dedup_type.value != "discrete" or self.dedup_method.value != "priority":
            priorities = None
        elif not priorities:
            ui.notify("Enter at least one priority value.", type="warning")
            return
        config = {
            "condition_cols": list(self.dedup_conditions.value or []) or None,
            "data_type": self.dedup_type.value, "method": self.dedup_method.value, "priority": priorities,
            "aggregation": self.dedup_aggregation.value if self.dedup_type.value == "continuous" else "mean",
            "log_transform": self.dedup_transform.value if self.dedup_type.value == "continuous" else "None",
            "dropna_conditions": bool(self.dedup_drop_na.value),
        }

        def action(pipeline: DiptoxPipeline, progress: Callable, status: Callable) -> None:
            status("Grouping and deduplicating records")
            pipeline.config_deduplicator(**config)
            pipeline.dataset_deduplicate(progress_callback=progress)
        await self._run_pipeline_job(
            "Deduplication", action,
            lambda pipeline: f"Deduplication complete: {len(pipeline.df):,} records; {len(pipeline.excluded_df):,} excluded.",
        )

    async def _run_substructure_search(self) -> None:
        if not self._require_data():
            return
        pattern = (self.search_pattern.value or "").strip()
        if not pattern:
            ui.notify("Enter a SMARTS or SMILES pattern.", type="warning")
            return
        is_smarts = bool(self.search_smarts.value)

        def action(pipeline: DiptoxPipeline, _progress: Callable, status: Callable) -> None:
            status("Searching substructures")
            pipeline.substructure_search(pattern, is_smarts=is_smarts)
        await self._run_pipeline_job("Substructure search", action, f"Search complete; added Substructure_{pattern}.")

    async def _run_atom_filter(self) -> None:
        if not self._require_data():
            return
        kwargs: dict[str, int] = {}
        if self.use_heavy.value:
            kwargs.update(min_heavy_atoms=int(self.min_heavy.value or 0), max_heavy_atoms=int(self.max_heavy.value or 0))
        if self.use_total.value:
            kwargs.update(min_total_atoms=int(self.min_total.value or 0), max_total_atoms=int(self.max_total.value or 0))
        if not kwargs:
            ui.notify("Enable at least one atom-count range.", type="warning")
            return

        def action(pipeline: DiptoxPipeline, _progress: Callable, status: Callable) -> None:
            status("Filtering records")
            pipeline.filter_by_atom_count(**kwargs)
        await self._run_pipeline_job("Atom-count filter", action, lambda pipeline: f"Filter complete: {len(pipeline.df):,} records remain.")

    def _recommended_columns(self) -> list[str]:
        pipeline = self.pipeline
        if pipeline.df is None:
            return []
        preferred = [
            *pipeline.df.attrs.get('_diptox_merged_columns', []),
            *pipeline.df.attrs.get('_diptox_condition_columns', {}).keys(),
            *pipeline.df.attrs.get('_diptox_condition_columns', {}).values(),
            pipeline.id_col, pipeline.name_col, pipeline.cas_col, "Is Valid", "Standardization Status",
            pipeline.smiles_col, "Canonical SMILES", pipeline.target_col, pipeline.unit_col,
            *((pipeline._current_dedup_config or {}).get('condition_cols') or []),
            "Unit Conversion Status", "Deduplication Strategy", "Deduplication Record Count",
            "Deduplication Source Rows", "Deduplication Input Values", "Structure Changed",
            "Is Multi-Component", "Final Contains Metal", "InChI",
        ]
        return [column for column in dict.fromkeys(preferred) if column and column in pipeline.df.columns]

    def _select_recommended_columns(self) -> None:
        self._export_selection_mode = 'recommended'
        self._set_export_columns(self._recommended_columns())

    def _export_columns_changed(self, _event=None) -> None:
        if not self._export_refreshing:
            self._export_selection_mode = 'custom'

    def _set_export_columns(self, values, options=None) -> None:
        self._export_refreshing = True
        try:
            if options is None:
                self.export_columns.set_value(values)
            else:
                self.export_columns.set_options(options, value=values)
        finally:
            self._export_refreshing = False

    def _clear_export_columns(self) -> None:
        self._export_selection_mode = 'custom'
        self._set_export_columns([])

    def _select_all_columns(self) -> None:
        columns = list(map(str, self.pipeline.df.columns)) if self.pipeline.df is not None else []
        self._export_selection_mode = 'all'
        self._set_export_columns(columns)

    @staticmethod
    def _temporary_export_path(stem: str, extension: str) -> Path:
        directory = Path(tempfile.gettempdir()) / "diptox-exports"
        directory.mkdir(parents=True, exist_ok=True)
        return directory / f"{stem}-{uuid.uuid4().hex[:8]}.{extension}"

    async def _export_results(self) -> None:
        if not self._require_data():
            return
        columns = list(self.export_columns.value or [])
        if not columns:
            ui.notify("Select at least one column.", type="warning")
            return
        extension = self.export_format.value or "csv"
        try:
            filename = self._download_filename(extension)
        except ValueError as error:
            ui.notify(str(error), type="warning")
            return
        output = self._temporary_export_path("diptox-processed", extension)
        source = self.pipeline

        def operation(_progress: Callable, status: Callable) -> Path:
            status("Writing result file")
            snapshot = _clone_pipeline(source)
            snapshot.save_results(str(output), columns=columns)
            return output
        result = await self._run_read_job("Exporting results", operation)
        if result:
            ui.download(result, filename=filename)
            ui.notify("Your download is ready.", type="positive")

    async def _export_excluded(self) -> None:
        if self.pipeline.excluded_df.empty:
            ui.notify("No excluded records are available.", type="warning")
            return
        try:
            filename = self._download_filename('csv', excluded=True)
        except ValueError as error:
            ui.notify(str(error), type="warning")
            return
        output = self._temporary_export_path("diptox-excluded", "csv")
        source = self.pipeline

        def operation(_progress: Callable, status: Callable) -> Path:
            status("Writing excluded records")
            snapshot = _clone_pipeline(source)
            snapshot.save_excluded_results(str(output))
            return output
        result = await self._run_read_job("Exporting excluded records", operation)
        if result:
            ui.download(result, filename=filename)
            ui.notify("Excluded records are ready.", type="positive")

    def _download_filename(self, extension: str, excluded: bool = False) -> str:
        name = (self.export_filename.value or '').strip()
        if Path(name).suffix.lower() in {'.csv', '.xlsx', '.txt', '.sdf', '.smi'}:
            name = name[: -len(Path(name).suffix)]
        if not name or name in {'.', '..'} or name.endswith('.') or re.search(r'[<>:"/\\|?*\x00-\x1f]', name):
            raise ValueError('Enter a file name without path separators or invalid filename characters.')
        return f"{name}{'-excluded' if excluded else ''}.{extension}"

    async def _undo(self) -> None:
        if not self.pipeline._history:
            ui.notify("No previous processing step is available.", type="warning")
            return

        def action(pipeline: DiptoxPipeline, _progress: Callable, status: Callable) -> None:
            status("Restoring previous step")
            if not pipeline.undo():
                raise ValueError("No previous processing step is available.")
        await self._run_pipeline_job("Undo", action, "Previous pipeline state restored.")

    def _refresh_data_controls(self) -> None:
        pipeline = self.pipeline
        frame = pipeline.df
        if frame is None:
            columns: list[str] = []
            self.header_dataset.set_text("Not loaded")
            self.sidebar_summary.set_text("No dataset loaded")
        else:
            columns = list(map(str, frame.columns))
            self.header_dataset.set_text(f"{len(frame):,} × {len(columns):,}")
            self.sidebar_summary.set_text(f"{len(frame):,} rows\n{len(columns):,} columns")
        pairs = frame.attrs.get('_diptox_condition_columns', {}) if frame is not None else {}
        merged_columns = frame.attrs.get('_diptox_merged_columns', []) if frame is not None else []
        condition_columns = [column for column in dict.fromkeys([*pipeline.source_columns, *pairs, *pairs.values(), *merged_columns]) if column in columns]
        conditions = [value for value in (self.dedup_conditions.value or []) if value in condition_columns]
        self.dedup_conditions.set_options(condition_columns, value=conditions)
        previous_export = [value for value in (self.export_columns.value or []) if value in columns]
        selected = (self._recommended_columns() if self._export_selection_mode == 'recommended'
                    else columns if self._export_selection_mode == 'all' else previous_export)
        self._set_export_columns(selected, options=columns)
        selected_column = self.value_filter_column.value if self.value_filter_column.value in columns else None
        self.value_filter_column.set_options(columns, value=selected_column)
        self._value_filter_column_changed()
        selected_merge = self.value_merge_column.value if self.value_merge_column.value in columns else None
        self.value_merge_column.set_options(columns, value=selected_merge)
        self._value_merge_column_changed()
        excluded_count = len(pipeline.excluded_df)
        self.excluded_summary.set_text(
            f"{excluded_count:,} excluded records are ready for download." if excluded_count else "No excluded records are available."
        )
        self.export_excluded_button.set_enabled(excluded_count > 0 and not self.store.busy)

        condition = self.unit_scope.value == 'Condition column'
        for control in (self.condition_value, self.condition_unit):
            control.set_options(columns, value=control.value if control.value in columns else None)
            control.set_enabled(condition)
        self.condition_log.set_enabled(condition)
        valid_units = bool(frame is not None and (condition or (pipeline.target_col in columns and pipeline.unit_col in columns)))
        units = self._detected_units()
        standard = (_optional_text(self.unit_standard.value) or '') if condition else (_optional_text(self.unit_standard.value) or (units[0] if units else ""))
        self.unit_standard.set_autocomplete(units)
        self.unit_standard.set_value(standard)
        mw_value = self.mw_column.value if self.mw_column.value in columns else None
        self.mw_column.set_options(columns, value=mw_value)
        self.unit_run_button.set_enabled(valid_units and not self.store.busy)
        if condition:
            self.unit_context.set_text('Condition column: select value/unit columns below; results are available for grouping and export.')
        elif valid_units:
            self.unit_context.set_text(
                f"Target: {pipeline.target_col}  ·  Unit column: {pipeline.unit_col}  ·  Detected: {', '.join(units) or 'none'}"
            )
        else:
            self.unit_context.set_text("Load data with target and unit columns to configure this step.")
        self.unit_rules_refresh.refresh()
        self.overview_refresh.refresh()
        self.undo_button.set_enabled(not self.store.busy and bool(pipeline._history))


_page_registered = False


def create_web_app() -> None:
    """Register the DiPTox page once for tests and entry points."""
    global _page_registered
    if _page_registered:
        return

    @ui.page("/")
    def index() -> None:
        DiptoxWebApp()

    _page_registered = True


def main(argv: Optional[list[str]] = None) -> int:
    """Start the local browser interface."""
    del argv
    multiprocessing.freeze_support()
    create_web_app()
    ui.run(
        host="127.0.0.1", port=int(os.environ.get("DIPTOX_PORT", "8080")),
        title=f"DiPTox · {FULL_NAME}", language="en-US", show=True,
        native=False, reload=False, show_welcome_message=False,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
