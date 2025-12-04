# diptox/__init__.py
from .core import DiptoxPipeline
from .chem_processor import ChemistryProcessor
from .web_request import WebService
from .data_io import DataHandler
from .data_deduplicator import DataDeduplicator
from .logger import LogManager
from .substructure_search import SubstructureSearcher
from .unit_processor import UnitProcessor

__all__ = ["DiptoxPipeline",
           "ChemistryProcessor",
           "WebService",
           "DataHandler",
           "DataDeduplicator",
           "LogManager",
           "SubstructureSearcher",
           "UnitProcessor"
           ]
__version__ = "1.3.1"
