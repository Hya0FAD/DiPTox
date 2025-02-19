# molppc/web_request.py
import requests
from typing import Optional, Tuple, List, Dict, Any, Set, Callable
import time
import re
from threading import Lock
from tqdm import tqdm
from rdkit import Chem
from concurrent.futures import ThreadPoolExecutor, as_completed
from .logger import log_manager
logger = log_manager.get_logger(__name__)


class DataSource:
    PUBCHEM = 'pubchem'
    CHEMSPIDER = 'chemspider'
    COMPTOX = 'comptox'
    CACTUS = 'cactus'


class WebService:
    """Handles all network request operations."""

    def __init__(self, source: str = DataSource.PUBCHEM,
                 interval: int = 0.3,
                 retries: int = 3,
                 delay: int = 30,
                 max_workers: int = 4,
                 batch_limit: int = 1500,
                 rest_duration: int = 300,
                 chemspider_api_key: Optional[str] = None,
                 comptox_api_key: Optional[str] = None):
        """
        Initializes the WebService class
        :param source: Data source interface
        :param interval: Time interval in seconds
        :param retries: Number of retry attempts on failure.
        :param delay: Delay between retries (in seconds).
        :param max_workers: Maximum number of concurrent requests.
        :param batch_limit: Number of requests before taking a break.
        :param rest_duration: Duration of the break in seconds.
        :param chemspider_api_key: Chemspider API key.
        :param comptox_api_key: Comptox API key.
        """
        self.source = source.lower()
        self.interval = interval
        self.retries = retries
        self.delay = delay
        self.max_workers = max_workers
        self.batch_limit = batch_limit
        self.rest_duration = rest_duration
        self.request_count = 0
        self._counter_lock = Lock()
        self._chemspider_api_key = chemspider_api_key
        self._comptox_api_key = comptox_api_key
        self._cas_pattern = re.compile(r'\d+-\d+-\d')
        self._valid_sources = {
            DataSource.PUBCHEM: self._fetch_via_pubchem,
            DataSource.CHEMSPIDER: self._fetch_via_chemspider,
            DataSource.COMPTOX: self._fetch_via_comptox,
            DataSource.CACTUS: self._fetch_via_cactus,
        }

    def get_properties(self,
                       identifiers: List[str],
                       properties: Set[str],
                       identifier_type: str = 'smiles') -> List[Dict[str, Optional[str]]]:
        logger.debug(f'Starting {self.source} query for {len(identifiers)} items')
        return self._batch_process(
            identifiers,
            lambda x: self._fetch_properties_wrapper(x, properties, identifier_type)
        )

    def _fetch_properties_wrapper(self,
                                  identifier: str,
                                  properties: Set[str],
                                  identifier_type: str) -> Dict[str, Optional[str]]:
        """Wrap the attribute acquisition logic to maintain a uniform retry mechanism"""
        fetch_fn = self._valid_sources.get(self.source)
        if not fetch_fn:
            raise ValueError(f"Invalid data source: {self.source}")
        try:
            result = fetch_fn(identifier, properties, identifier_type)
            time.sleep(self.interval)
        except Exception as e:
            raise e

        return {prop: result.get(prop) for prop in properties}

    def _batch_process(self, smiles_list: List[str],
                       fetch_func: Callable[[str], Dict]) -> List[Any]:
        """General framework for batch processing."""
        logger.info(f"Initializing batch processing with {self.max_workers} workers")

        results = [None] * len(smiles_list)
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(self._retry_wrapper, fetch_func, smiles, index): index
                       for index, smiles in enumerate(smiles_list)}

            for future in tqdm(as_completed(futures), total=len(futures)):
                index, result = future.result()
                results[index] = result

        logger.info(f"Completed batch processing with {len(smiles_list)} items")
        return results

    def _retry_wrapper(self, func: callable,
                       identifier: str,
                       index: int) -> Tuple[int, Any]:
        """Wrapper that implements retry logic for failed requests."""
        for attempt in range(1, self.retries + 1):
            try:
                result = func(identifier)
                return index, result
            except Exception:
                logger.warning(f"Attempt {attempt}/{self.retries} failed for {identifier}")
                if attempt < self.retries:
                    time.sleep(self.delay * (2 ** (attempt - 1)))
        logger.error(f"Permanent failure for {identifier}")
        return index, {prop: None for prop in getattr(func, 'properties', set())}

    def _increment_request_count(self):
        """
        Increment the request count and check if a break is needed.
        If the number of requests exceeds the batch limit, pause the execution for a specified duration.
        """
        with self._counter_lock:
            self.request_count += 1
            if self.request_count >= self.batch_limit:
                logger.warning(f"Reached request limit {self.batch_limit}, entering cooldown")
                time.sleep(self.rest_duration)
                self.request_count -= self.batch_limit

    def _fetch_via_pubchem(self, identifier: str,
                           properties: Set[str],
                           identifier_type: str) -> Dict[str, Optional[str]]:
        """Fetch through PubChem's API"""
        result = {}
        try:
            from pubchempy import get_compounds, Compound
            if identifier_type == 'smiles':
                compounds = get_compounds(identifier, 'smiles')
            elif identifier_type == 'cas':
                compounds = get_compounds(identifier, 'name')
            else:
                logger.error(f"Invalid identifier type {identifier_type}")
                raise

            # if isinstance(compounds, dict) and 'Fault' in compounds:
            #     if 'ServerBusy' in compounds.get('Fault', ''):
            #         raise requests.exceptions.HTTPError()

            self._increment_request_count()
            if not compounds:
                return {prop: None for prop in properties}

            compound = compounds[0]
            result['cid'] = str(compound.cid)

            if 'cas' in properties:
                result['cas'] = self._parse_pubchem_cas_direct(compound)
            if 'iupac' in properties:
                result['iupac'] = compound.iupac_name
            if 'smiles' in properties:
                result['smiles'] = compound.canonical_smiles

            return result

        except Exception as e:
            logger.error(f"PubChem query failed for {identifier}: {str(e)}")
            raise

    def _parse_pubchem_cas_direct(self, compound) -> Optional[str]:
        """Resolve the CAS number directly from the Compound object"""
        try:
            for synonym in compound.synonyms:
                if self._cas_pattern.match(synonym):
                    return synonym
            return None
        except Exception as e:
            logger.error(f"CAS parsing failed: {str(e)}")
            return None

    def _fetch_via_chemspider(self, identifier: str,
                              properties: Set[str],
                              identifier_type: str) -> Dict[str, Optional[str]]:
        """Fetch through ChemSpider's API"""
        if not self._chemspider_api_key:
            logger.error("ChemSpider API key required")
            raise

        result = {}
        try:
            from chemspipy import ChemSpider
            cs = ChemSpider(self._chemspider_api_key)
            search_results = cs.search(identifier)
            self._increment_request_count()

            for _ in range(0, self.retries):
                if search_results.status in {'Complete', 'Failed'}:
                    break
                else:
                    time.sleep(2)
                    continue
            if search_results.status != 'Complete':
                return {prop: None for prop in properties}

            compound = search_results[0]
            print(compound)
            if 'cas' in properties:
                result['cas'] = compound.get('cas')
            if 'iupac' in properties:
                result['iupac'] = compound.get('iupac_name')
            if 'smiles' in properties:
                result['smiles'] = compound.smiles

            return result

        except Exception as e:
            logger.error(f"ChemSpider query failed: {str(e)}")
            raise

    def _fetch_via_comptox(self, identifier: str,
                           properties: Set[str],
                           identifier_type: str) -> Dict[str, Optional[str]]:
        """Fetch through Comptox's API"""
        if not self._comptox_api_key:
            logger.error("CompTox API key required")
            raise

        result = {}
        try:
            import ctxpy as ctx
            comptox = ctx.Chemical(x_api_key=self._comptox_api_key)
            search_results = comptox.search(by='equals', word=identifier)
            self._increment_request_count()

            if not search_results:
                return {prop: None for prop in properties}

            dtxsid = search_results[0].get('dtxsid')
            detail = comptox.details(by='dtxsid', word=dtxsid)
            self._increment_request_count()
            if detail.success and detail.records:
                record = detail.records[0]
                if 'cas' in properties:
                    result['cas'] = record.get('casrn')
                if 'iupac' in properties:
                    result['iupac'] = record.get('iupac_name')
                if 'smiles' in properties:
                    result['smiles'] = record.get('smiles')

            return result

        except Exception as e:
            logger.error(f"CompTox query failed: {str(e)}")
            raise

    def _fetch_via_cactus(self, identifier: str,
                          properties: Set[str],
                          identifier_type: str) -> Dict[str, Optional[str]]:
        """Fetch through Cactus's API"""
        result = {}
        try:
            if identifier_type == 'smiles':
                mol = Chem.MolFromSmiles(identifier)
                if not mol:
                    logger.error("Invalid SMILES")
                    raise
                identifier = Chem.MolToInchiKey(mol)

            base_url = f"https://cactus.nci.nih.gov/chemical/structure/{identifier}"

            if 'cas' in properties:
                response = requests.get(f"{base_url}/cas", timeout=10)
                self._increment_request_count()
                result['cas'] = response.text.strip() if response.ok else None
            if 'iupac' in properties:
                response = requests.get(f"{base_url}/iupac_name", timeout=10)
                self._increment_request_count()
                result['iupac'] = response.text.strip() if response.ok else None
            if 'smiles' in properties:
                response = requests.get(f"{base_url}/smiles", timeout=10)
                self._increment_request_count()
                result['smiles'] = response.text.strip() if response.ok else None

            return result

        except Exception as e:
            logger.error(f"Cactus query failed: {str(e)}")
            raise
