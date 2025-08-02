# Arthor

__all__ = ['QueryArthor']

import json
import logging
import os
import warnings
from datetime import datetime
from pathlib import Path
from time import sleep
from typing import List, Optional, Dict, Union

import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry


class QueryArthor:
    """
    Query class for Arthorian Quest.
    It queries arthor.docking.org via ``.retrieve`` or ``.batch_retrieve``

    See https://arthor.docking.org/api.html for the API endpoints used

    .. code-block:: python
        Query().retrieve('[CH3]-[CH2X4]', ['BB-50-22Q1'])

        # For batch queries:
        Query().batch_retrieve(['[CH3]-[CH2X4]', '[OH]-[CH2]'], ['BB-50-22Q1'])
    """

    enamine_dbs = ['BB-ForSale-22Q1', 'MADE-BB-23Q1-770M', 'REAL-Database-22Q1']

    def __init__(self, base_url: str = 'https://arthor.docking.org/',
                 cache_dir: Optional[Union[str, Path]] = None,
                 max_retries: int = 3,
                 backoff_factor: float = 1.0):
        self.base_url = base_url
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(exist_ok=True)

        # Configure session with retry strategy
        self.session = requests.Session()
        retry_strategy = Retry(
            total=max_retries,
            backoff_factor=backoff_factor,
            status_forcelist=[500, 502, 503, 504],
            allowed_methods=["GET"]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

    @property
    def dbs(self):
        return pd.DataFrame(self.session.get(self.base_url + 'dt/data').json())

    def create_empty_metadata_row(self, query: str, query_inchi: str, data: dict, length: int,
                                  dbname: str) -> pd.DataFrame:
        """Creates a single row DataFrame with metadata when no matches are found"""
        try:
            arthor_source = data['arthor.source']
        except KeyError:
            arthor_source = dbname
        return pd.DataFrame([{
            'arthor.rank': None,
            'arthor.index': None,
            'smiles': None,
            'identifier': None,
            'arthor.source': arthor_source,
            'recordsTotal': data['recordsTotal'],
            'recordsFiltered': data['recordsFiltered'],
            'hasMore': data['hasMore'],
            'query': query,
            'query_inchi': query_inchi,
            'query_length': 0,
            'N_RB': None,
            'N_HA': None,
            'mol': None
        }])

    def retrieve(self, query: str, dbnames: List[str], search_type='Substructure', length=1_000):
        """
        Returns a dataframe of the results of the query,
        with fields:

        * N_RB: number of rotatable bonds
        * N_HA: number of heavy atoms

        :param query: SMARTS query
        :param dbnames: list of names (see self.dbs)
        :return:
        """
        try:
            dbname: str = ','.join(dbnames)
            query_inchi: Optional[str] = None

            # Validate and process query
            if isinstance(query, Chem.Mol):
                query = Chem.MolToSmarts(query)
                query_inchi = Chem.MolToInchiKey(query)
            if isinstance(query, str) and query_inchi is None:
                mol = Chem.MolFromSmiles(query)
                if mol is None:
                    raise ValueError(f"Invalid SMILES query: {query}")
                query_inchi = Chem.MolToInchiKey(mol)

            # Use session with retry logic and timeout
            response = self.session.get(
                self.base_url + f'/dt/{dbname}/search',
                params={
                    'query': query,
                    'type': search_type,
                    'length': length
                },
                timeout=30
            )

            if response.status_code == 503:
                raise ConnectionError('Arthor unavailable. cf. https://arthor.docking.org/')

            response.raise_for_status()
            data: dict = response.json()

            # Handle empty or invalid responses
            if data.get("message", '') == "SMARTS query is always false!":
                warnings.warn(f"SMARTS query {query} is always false")
                return self.create_empty_metadata_row(query=query,
                                                      query_inchi=query_inchi,
                                                      data=data,
                                                      length=length,
                                                      dbname=dbname)
            if data.get('warning', ''):
                warnings.warn(data['warning'])
            if not data.get('recordsTotal', False):
                warnings.warn(f"SMARTS query {query} returned no matches")
                return self.create_empty_metadata_row(query=query,
                                                      query_inchi=query_inchi,
                                                      data=data,
                                                      length=length,
                                                      dbname=dbname)

            matches = pd.DataFrame(data['data'],
                                   columns=['arthor.rank', 'arthor.index', 'smiles', 'identifier', 'arthor.source'])

            if len(matches) == 0:  # empty
                return self.create_empty_metadata_row(query=query,
                                                      query_inchi=query_inchi,
                                                      data=data,
                                                      length=length,
                                                      dbname=dbname)

            matches['arthor.source'] = matches['arthor.source'].apply(lambda x: x.replace('\t', ''))
            # add metadata from query
            matches['recordsTotal'] = data['recordsTotal']
            matches['recordsFiltered'] = data['recordsFiltered']
            matches['hasMore'] = data['hasMore']
            matches['query'] = query
            matches['query_inchi'] = query_inchi
            matches['query_length'] = length

            # Remove duplicates before heavy computations
            matches = matches.drop_duplicates('arthor.index')

            # Add molecule information with error handling
            PandasTools.AddMoleculeColumnToFrame(matches, 'smiles', 'mol', includeFingerprints=True)
            matches = matches.loc[~matches.mol.isnull()]

            # Calculate properties
            matches['N_RB'] = matches.mol.apply(AllChem.CalcNumRotatableBonds)
            matches['N_HA'] = matches.mol.apply(AllChem.CalcNumHeavyAtoms)

            return matches.sort_values('N_HA').reset_index(drop=True)

        except requests.exceptions.RequestException as e:
            logging.error(f"Network error during query {query}: {str(e)}")
            raise
        except Exception as e:
            logging.error(f"Error processing query {query}: {str(e)}")
            raise

    def _is_batch_complete(self, batch_id: str, batch_file: str, query_inchis: List[str],
                           completed_batches: Dict) -> bool:
        """Helper method to check if a batch is complete"""
        if batch_id in completed_batches and completed_batches[batch_id]:
            if os.path.exists(batch_file):
                try:
                    batch_df = pd.read_pickle(batch_file)
                    return all(inchi in batch_df.query_inchi.values for inchi in query_inchis)
                except Exception as e:
                    logging.warning(f"Could not read batch file {batch_file}: {str(e)}")
        return False

    def _update_progress(self, progress_file: str, batch_id: str,
                         completed_batches: Dict, status: bool) -> None:
        """Helper method to update progress safely"""
        try:
            completed_batches[batch_id] = status
            with open(progress_file, 'w') as f:
                json.dump(completed_batches, f)
        except Exception as e:
            logging.error(f"Failed to update progress file: {str(e)}")

    def batch_retrieve(self,
                       queries: List[str],
                       dbnames: List[str],
                       search_type: str = 'Substructure',
                       length: int = 10_000,
                       sleep_time: float = 5.0,
                       continue_on_error: bool = True,
                       batch_size: int = 10) -> str:
        """
        Perform batch retrieval of multiple queries with progress tracking and error handling.

        Parameters
        ----------
        queries : List[str]
            List of SMARTS or SMILES queries
        dbnames : List[str]
            List of database names to search
        search_type : str, optional
            Type of search ('SMARTS' or 'Substructure' or 'Similarity'), by default 'Substructure'
        length : int, optional
            Maximum number of results per database per query, by default 10_000.
        sleep_time : float, optional
            Time to wait between queries in seconds, by default 5.0
        continue_on_error : bool, optional
            Whether to continue processing on error, by default True
        batch_size : int, optional
            Number of queries to save in intermittent dataframes, by default 10

        Returns
        -------
        str
            Path to the cache directory where batched dataframes are saved
        """
        # Setup logging
        log_file = os.path.join(
            self.cache_dir if self.cache_dir else Path.cwd(),
            f"arthorian_quest_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        )
        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )

        if not self.cache_dir:
            self.cache_dir = Path.cwd()
            logging.warning(f"No cache directory set. Using current directory: {self.cache_dir}")

        # Initialize progress tracking
        progress_file = os.path.join(self.cache_dir, "batch_progress.json")
        completed_batches = {}
        if os.path.exists(progress_file):
            try:
                with open(progress_file, 'r') as f:
                    completed_batches = json.load(f)
            except json.JSONDecodeError:
                logging.warning("Corrupted progress file. Starting fresh.")

        # Input validation
        if not all(isinstance(q, str) for q in queries):
            logging.error("Queries must be a list of strings!")
            logging.error(f"Queries provided: {queries}")
            raise TypeError

        if not all(isinstance(db, str) for db in dbnames):
            logging.error("dbnames must be a list of strings!")
            logging.error(f"dbnames provided: {dbnames}")
            raise TypeError

        # Create batches of queries
        query_batches = [queries[i:i + batch_size] for i in range(0, len(queries), batch_size)]

        for batch_idx, query_batch in enumerate(tqdm(query_batches, desc="Processing batches")):
            batch_id = f"batch_{batch_idx}"
            batch_file = os.path.join(self.cache_dir, f"batch_{batch_idx}.pkl.gz")

            try:
                # Convert queries to InChI keys with error handling
                query_inchis = []
                for q in query_batch:
                    mol = Chem.MolFromSmiles(q)
                    if mol is None:
                        raise ValueError(f"Invalid SMILES in batch {batch_idx}: {q}")
                    query_inchis.append(Chem.MolToInchiKey(mol))

                # Skip if batch is already completed
                if self._is_batch_complete(batch_id, batch_file, query_inchis, completed_batches):
                    logging.info(f"Skipping completed batch {batch_idx}")
                    continue

                batch_results = []
                for query in query_batch:
                    for dbname in dbnames:
                        try:
                            df = self.retrieve(query, [dbname], search_type, length)
                            if df is not None and not df.empty:
                                batch_results.append(df)
                                logging.info(f"Successfully processed query: {query} {dbname} ({len(df)} results)")
                        except Exception as e:
                            error_msg = f"Error processing query {query} {dbname}: {str(e)}"
                            logging.error(error_msg)
                            if not continue_on_error:
                                raise
                    sleep(sleep_time)

                if batch_results:
                    batch_df = pd.concat(batch_results, ignore_index=True)
                    batch_df.to_pickle(batch_file)
                    self._update_progress(progress_file, batch_id, completed_batches, True)
                    logging.info(f"Completed and saved batch {batch_idx}")

            except Exception as e:
                error_msg = f"Error processing batch {batch_idx}: {str(e)}"
                logging.error(error_msg)
                self._update_progress(progress_file, batch_id, completed_batches, False)
                if not continue_on_error:
                    raise Exception(error_msg)

        return self.cache_dir

    def get_batch_statistics(self) -> Dict:
        """
        Get statistics about batch processing progress.
        Only available if cache_dir was set.

        Returns
        -------
        Dict
            Dictionary containing progress statistics
        """
        if not self.cache_dir:
            return {"error": "No cache directory set"}

        progress_file = self.cache_dir / "batch_progress.json"
        if not progress_file.exists():
            return {"error": "No batch processing history found"}

        with open(progress_file, 'r') as f:
            completed = json.load(f)

        return {
            "total_processed": len(completed),
            "successful": sum(1 for v in completed.values() if v),
            "failed": sum(1 for v in completed.values() if not v)
        }
