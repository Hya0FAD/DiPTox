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
import logging
from .logger import log_manager
logger = log_manager.get_logger(__name__)


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
    _CTXPY_AVAILABLE = True
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
                 force_api_mode: bool = False):
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
        self._chemspider_api_key = chemspider_api_key
        self._comptox_api_key = comptox_api_key
        self._cas_api_key = cas_api_key

        self.session = requests.Session()
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
                logger.info(f"The request limit ({self.batch_limit}) has been reached. Pausing for {self.rest_duration} seconds...")
                time.sleep(self.rest_duration)
                self.request_count = 0

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
        default_error_result = {prop: None for prop in properties}
        default_error_result['Data_Source'] = 'Error'

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
                    logger.error(f"An error occurred while processing the identifier '{identifiers}': {e}")
                    results_map[identifier] = default_error_result
                completed_count += 1
                if progress_callback:
                    progress_callback(completed_count, total_tasks)
        final_results = []
        for original_id in identifiers:
            if is_valid_id(original_id) and original_id in results_map:
                final_results.append(results_map[original_id])
            else:
                empty_res = {prop: None for prop in properties}
                if is_valid_id(original_id):
                    empty_res['Data_Source'] = 'Error'
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
        contributing_sources = set()

        for source in self.sources:
            if not needed_props:
                break
            fetch_fn = self._fetch_functions.get(source)
            if not fetch_fn:
                continue

            try:
                source_result = self._retry_wrapper(fetch_fn, identifier, needed_props, identifier_type)
                found_this_round = set()
                for prop, value in source_result.items():
                    if prop in needed_props and value is not None:
                        final_result[prop] = value
                        found_this_round.add(prop)
                if found_this_round:
                    contributing_sources.add(source)
                    needed_props -= found_this_round

            except Exception as e:
                logger.warning(f"Data source '{source}' failed to query '{identifier}': {e}. Trying the next one.")

        if contributing_sources:
            final_result['Data_Source'] = ", ".join(sorted(list(contributing_sources)))

        return final_result

    def _retry_wrapper(self,
                       func: Callable,
                       identifier: str,
                       properties: Set[str],
                       identifier_type: str) -> Dict[str, Any]:
        """Wrapper that implements retry logic for failed requests."""
        for attempt in range(self.retries):
            try:
                time.sleep(self.interval)
                return func(identifier, properties, identifier_type)
            except requests.exceptions.HTTPError as e:
                if e.response is not None and e.response.status_code:
                    status_code = e.response.status_code
                    if 400 <= status_code < 500:
                        return {prop: None for prop in properties}
                    elif 500 <= status_code < 600:
                        logger.warning(
                            f"The {attempt + 1}/{self.retries} attempt failed (server error {status_code}): {identifier} at {func.__name__}")
                else:
                    logger.warning(
                        f"On the {attempt + 1}/{self.retries} attempt, the operation failed (HTTP error, no status code): {identifier} at {func.__name__}: {e}")
            except requests.exceptions.RequestException as e:
                logger.warning(
                    f"On the {attempt + 1}/{self.retries} attempt, the operation failed (network error): {identifier} at {func.__name__}: {e}")
            if attempt < self.retries - 1:
                time.sleep(self.delay * (attempt + 1))
            else:
                logger.error(f"All {self.retries} attempts failed: {identifier} at {func.__name__}")
        return {prop: None for prop in properties}

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
        compounds = pubchempy.get_compounds(identifier, id_type_map[id_type])
        result = {}
        if compounds:
            c = compounds[0]
            if 'smiles' in properties:
                result['smiles'] = c.canonical_smiles
            if 'iupac' in properties:
                result['iupac'] = c.iupac_name
            if 'mw' in properties:
                result['mw'] = c.molecular_weight
            if c.synonyms is not None:
                if 'cas' in properties:
                    result['cas'] = next((syn for syn in c.synonyms if self._validate_and_clean_cas(syn)), None)
                if 'name' in properties and c.synonyms:
                    result['name'] = c.synonyms[0]
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

            try:
                prop_url = f"{base_url}/compound/{id_namespace}/{quote(id_value)}/property/{prop_str}/JSON"
                self._increment_request_count()
                response = self.session.get(prop_url, timeout=15)
                if response.status_code == 200:
                    data = response.json()
                    if props := data.get('PropertyTable', {}).get('Properties', []):
                        vals = props[0]
                        if 'smiles' in properties: parsed_result['smiles'] = vals.get('CanonicalSMILES') or vals.get(
                            'ConnectivitySMILES')
                        if 'iupac' in properties: parsed_result['iupac'] = vals.get('IUPACName')
                        if 'mw' in properties: parsed_result['mw'] = vals.get('MolecularWeight')
            except requests.exceptions.RequestException:
                logger.debug(f"Failed to obtain the single-valued attribute for '{id_value}'")
                return None

            if 'name' in properties or 'cas' in properties:
                try:
                    syn_url = f"{base_url}/compound/{id_namespace}/{quote(id_value)}/synonyms/JSON"
                    self._increment_request_count()
                    syn_response = self.session.get(syn_url, timeout=15)
                    if syn_response.status_code == 200:
                        syn_data = syn_response.json()
                        if synonyms := syn_data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym'):
                            if 'name' in properties:
                                parsed_result['name'] = synonyms[0]
                            if 'cas' in properties:
                                parsed_result['cas'] = next(
                                    (syn for syn in synonyms if self._validate_and_clean_cas(syn)), None)
                except requests.exceptions.RequestException:
                    logger.debug(f"Failed to obtain synonyms for '{id_value}'")

            return parsed_result if parsed_result else None

        result = get_and_parse(identifier, id_type_map[id_type])
        if result:
            return result

        try:
            cid_url = f"{base_url}/compound/{id_type_map[id_type]}/{quote(identifier)}/cids/JSON"
            self._increment_request_count()
            cid_response = self.session.get(cid_url, timeout=15)
            if cid_response.status_code == 200:
                if cids := cid_response.json().get('IdentifierList', {}).get('CID'):
                    result = get_and_parse(cids[0], 'cid')
                    if result:
                        return result
        except requests.exceptions.RequestException:
            logger.warning(f"Failed to query PubChem through CID fallback for '{identifier}'")

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
                return {}
            time.sleep(2)
        else:
            return {}

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

        for _ in range(3):
            self._increment_request_count()
            status_response = self.session.get(status_url, headers=headers_get, timeout=15)
            status_response.raise_for_status()
            status = status_response.json().get('status')

            if status == 'Complete':
                break
            if status == 'Failed':
                return {}
            time.sleep(2)
        else:
            logger.warning(f"ChemSpider search for '{identifier}' timed out.")
            return {}

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
                logger.warning(f"Unable to parse SMILES to generate InChIKey: {smiles}")
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
                logger.warning(f"Unable to parse SMILES to generate InChI: {smiles}")
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
        search_results = comptox.search(by='equals', word=search_identifier)

        if not search_results or not isinstance(search_results, list):
            return {}

        dtxsid = search_results[0].get('dtxsid')
        if not dtxsid:
            return {}

        self._increment_request_count()
        detail = comptox.details(by='dtxsid', word=dtxsid)

        result = {}
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
            props = mol.get('molecule_properties', {})
            if 'smiles' in properties:
                result['smiles'] = mol.get('molecule_structures', {}).get('canonical_smiles')
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
        if id_type == 'smiles':
            identifier = self._get_inchi_from_smiles(identifier)
            if not identifier:
                return {}
        elif id_type not in ['cas', 'name']:
            return {}
        base_url = "https://commonchemistry.cas.org/api"
        headers = {
            'Authorization': self._cas_api_key
        }
        search_url = f"{base_url}/search"
        params = {'q': identifier}

        self._increment_request_count()
        search_response = self.session.get(search_url, params=params, headers=headers, timeout=15)

        if search_response.status_code == 404:
            return {}
        search_response.raise_for_status()

        search_data = search_response.json()
        if not search_data.get('results'):
            return {}

        cas_rn = search_data['results'][0].get('rn')
        if not cas_rn:
            return {}

        detail_url = f"{base_url}/detail"
        params = {'cas_rn': cas_rn}

        self._increment_request_count()
        detail_response = self.session.get(detail_url, params=params, timeout=15)

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
            if synonyms := detail_data.get('synonyms'):
                result['name'] = synonyms[0]
            else:
                result['name'] = detail_data.get('name')

        return result