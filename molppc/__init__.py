# molppc/__init__.py
from .core import MolecularProcessor
from .chem_processor import ChemistryProcessor
from .web_request import WebService
from .data_io import DataHandler
from .data_deduplicator import DataDeduplicator
from .logger import LogManager

__all__ = ["MolecularProcessor", "ChemistryProcessor", "WebService", "DataHandler", "DataDeduplicator", "LogManager"]
__version__ = "0.3.1"
