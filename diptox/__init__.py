"""DiPTox's public API, loaded on demand for lightweight command discovery."""

from importlib import import_module

__all__ = ["DiptoxPipeline",
           "ChemistryProcessor",
           "WebService",
           "DataHandler",
           "DataDeduplicator",
           "LogManager",
           "SubstructureSearcher",
           "UnitProcessor"
           ]
__version__ = "1.1.0"

_PUBLIC_MODULES = {
    "DiptoxPipeline": "core",
    "ChemistryProcessor": "chem_processor",
    "WebService": "web_request",
    "DataHandler": "data_io",
    "DataDeduplicator": "data_deduplicator",
    "LogManager": "logger",
    "SubstructureSearcher": "substructure_search",
    "UnitProcessor": "unit_processor",
}


def __getattr__(name):
    module = _PUBLIC_MODULES.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f".{module}", __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
