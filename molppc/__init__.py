# molppc/__init__.py
from .core import MolecularProcessor
from .chem_processor import ChemistryProcessor
from .web_request import WebService
from .data_io import DataHandler
from .data_deduplicator import DataDeduplicator
from .logger import LogManager
from .substructure_search import SubstructureSearcher

__all__ = ["MolecularProcessor",
           "ChemistryProcessor",
           "WebService",
           "DataHandler",
           "DataDeduplicator",
           "LogManager",
           "SubstructureSearcher"
           ]
__version__ = "0.12.0"
