#!/usr/bin/env python3
"""
Arthor.py

This module contains the functionality for an Arthor search.
"""
import json
import logging
import os
import random
import sys
import time
from typing import (Any, List, Dict, Tuple, Optional)

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from rdkit import Chem

import syndirella.utils.fairy as fairy
from syndirella.database.DatabaseSearch import DatabaseSearch
from syndirella.utils.error import APIRetryLimitExceeded


class Arthor(DatabaseSearch):
    """
    This class contains information about the Arthor search. It will perform the Arthor search using the
    perform_database_search function. It will also store the results of the Arthor search as a .csv file.
    """

    def __init__(self):
        super().__init__()
        self.url = os.environ.get("ARTHOR_API_URL", "https://arthor.docking.org")
        self.api_key = os.environ.get("ARTHOR_API_KEY", None)  # Arthor may not require API key
        self.logger = logging.getLogger(f"{__name__}")
        
        # Configure session with retry strategy
        self.session = requests.Session()
        retry_strategy = requests.adapters.Retry(
            total=3,
            backoff_factor=1.0,
            status_forcelist=[500, 502, 503, 504],
            allowed_methods=["GET"]
        )
        adapter = requests.adapters.HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

    def perform_database_search(self,
                                reactant: Chem.Mol,
                                reaction_name: str,
                                search_type: str = "superstructure",
                                vendors: list[str] = ['enamine_real', 'stock', 'mcule', 'zinc']) -> List[str] | None:
        """
        This function is used to perform the Arthor search using the database_search_function.
        """
        # 1. Get additional similar reactant if reaction is one with additional reactants
        reactant_filters = fairy.load_reactant_filters()
        reactants: List[str] = fairy.find_similar_reactants(reactant=reactant,
                                                            reaction_name=reaction_name,
                                                            reactant_filters=reactant_filters)
        # 2. Perform the search for all
        hits_all: List[Tuple[str, Tuple[str, str] | None]] = []
        for smiles in reactants:
            if search_type == "superstructure":
                hits: List[Tuple[str, Tuple[str, str] | None]] | None = self.perform_superstructure_search(smiles, vendors=vendors)
                if hits is None:  # if the API query failed
                    return None
                self.logger.info(f'Found {len(hits)} hits for {smiles} before filtering.')
                hits_all.extend(hits)
        filtered_hits: List[str] = fairy.filter_molecules(hits=hits_all)
        return filtered_hits

    def perform_superstructure_search(self,
                                      smiles: str,
                                      queryThirdPartyServices: bool = False,
                                      keep_catalogue: bool = False,
                                      vendors: list[str] = ['all']) -> List[Tuple[str, Tuple[str, str] | None]] | None:
        """
        This function is used to perform the Arthor superstructure search.
        """
        self.logger.info(f"Running superstructure search for {smiles}.")
        self.logger.info(f"Querying third party services: {queryThirdPartyServices}")
        self.logger.info(f"Vendor list: {vendors}")
        if not isinstance(smiles, str):
            self.logger.error("Smiles must be a string.")
            raise TypeError("Smiles must be a string.")
        
        # Convert vendors list to Arthor database names
        arthor_dbs = self._convert_vendors_to_arthor_dbs(vendors)
        
        # Arthor uses GET requests with query parameters, not POST
        hits_all = []
        for dbname in arthor_dbs:
            try:
                self.logger.info(f"Querying Arthor database {dbname} for {smiles}.")
                response = self.session.get(
                    f'{self.url}/dt/{dbname}/search',
                    params={
                        'query': smiles,
                        'type': 'Substructure',
                        'length': 1000
                    },
                    timeout=30
                )
                
                if response.status_code == 503:
                    self.logger.warning(f'Arthor unavailable for database {dbname}')
                    continue
                    
                response.raise_for_status()
                data = response.json()
                
                if data.get('data'):
                    # data['data'] is a list of lists, and data['header'] is the list of column names
                    rows = data['data']
                    header = data.get('header')
                    if header:
                        self.logger.debug(f"Header: {header}")
                        # Store the raw data with header for processing
                        hits_all.append({
                            'data': rows,
                            'header': header,
                            'database': dbname
                        })
                    else:
                        self.logger.warning(f"No header found in response from {dbname}")

            except Exception as e:
                self.logger.error(f"Error querying database {dbname}: {str(e)}")
                continue
        
        hits_info: List[Tuple[str, Tuple[str, str] | None]] | None = self.structure_output(hits_all,
                                                                                           query_smiles=smiles,
                                                                                           keep_catalogue=keep_catalogue)

        return hits_info

    def _convert_vendors_to_arthor_dbs(self, vendors: list[str]) -> list[str]:
        """
        Convert vendor names to Arthor database names.
        """
        vendor_to_db_mapping = {
            # Enamine databases (actual available databases)
            'enamine_real': 'REAL-Database-22Q1',

            # Stock databases
            'stock': 'In-Stock-19Q4-14.1M',
            'chemspace': 'ChemSpace-SC-Stock-Mar2022-346K',
            
            # Mcule databases
            'mcule': 'Mcule-22Q1-8.7M',
            'mcule_bb': 'Mcule-BB-22Q1-2.1M',
            'mcule_full': 'Mcule-Full-22Q1-60M',
            'mcule_v': 'Mcule-V-22Q1-51M',
            'mcule_ultimate': 'Mcule-Ultimate-20Q2-126M',
            'mcule_purchasable': 'mcule_purchasable_virtual_230121',
            
            # ZINC databases
            'zinc_all': 'ZINC-All-19Q4-1.4B',
            'zinc_interesting': 'ZINC-Interesting-19Q4-307K',
            'zinc_on_demand': 'ZINC-On-Demand-19Q4-311M',
            'zinc_for_sale': 'ZINC20-ForSale-22Q1',
            
            # Grouped vendor databases (all major vendors)
            'enamine': 'REAL-Database-22Q1',
            'mcule': 'Mcule-22Q1-8.7M,Mcule-BB-22Q1-2.1M,Mcule-Full-22Q1-60M,Mcule-V-22Q1-51M,Mcule-Ultimate-20Q2-126M,mcule_purchasable_virtual_230121',
            'zinc': 'ZINC-All-19Q4-1.4B,ZINC-Interesting-19Q4-307K,ZINC-On-Demand-19Q4-311M,ZINC20-ForSale-22Q1',
            
            # Default mappings for backward compatibility
            'all': 'ZINC20-ForSale-22Q1,REAL-Database-22Q1,In-Stock-19Q4-14.1M,Mcule-22Q1-8.7M,Mcule-BB-22Q1-2.1M,ZINC-All-19Q4-1.4B,ZINC-Interesting-19Q4-307K,ZINC-On-Demand-19Q4-311M'
        }
        
        arthor_dbs = []
        for vendor in vendors:
            if vendor in vendor_to_db_mapping:
                db_names = vendor_to_db_mapping[vendor].split(',')
                for db_name in db_names:
                    if db_name not in arthor_dbs:
                        arthor_dbs.append(db_name)
        
        if not arthor_dbs:
            # Default to commonly available chemicals (In-Stock and ZINC20)
            arthor_dbs = ['In-Stock-19Q4-14.1M', 'ZINC20-ForSale-22Q1']
        
        return arthor_dbs
    
    def get_available_databases(self) -> list[str]:
        """
        Get the list of available databases from Arthor API.
        
        Returns
        -------
        list[str]
            List of available database names
        """
        try:
            response = self.session.get(f'{self.url}/dt/data')
            response.raise_for_status()
            databases = response.json()
            return [db['displayName'] for db in databases]
        except Exception as e:
            self.logger.error(f"Error fetching available databases: {str(e)}")
            return []
    


    def structure_output(self, 
                         hits: List[Dict] | None, 
                         query_smiles: str, 
                         keep_catalogue: bool = False) -> List[Tuple[str, Tuple[str, str] | None]] | None:
        """
        Formats output into a list of tuples with smiles and catalogue.
        """
        hits_info: List[Tuple[str, Tuple[str, str] | None]] = []
        if hits is None:
            self.logger.critical(
                f"Error with API output, returning None for superstructure search for {query_smiles}!")
            return None
        
        for hit_batch in hits:
            if not isinstance(hit_batch, dict) or 'data' not in hit_batch or 'header' not in hit_batch:
                self.logger.warning(f"Invalid hit batch structure: {hit_batch}")
                continue
                
            rows = hit_batch['data']
            header = hit_batch['header']
            database = hit_batch.get('database', 'unknown')
            
            # Find the column indices for the required fields
            try:
                smiles_idx = header.index('SMILES')
                identifier_idx = header.index('Identifier')
                source_idx = header.index('arthor.source') if 'arthor.source' in header else None
            except ValueError as e:
                self.logger.warning(f"Required column not found in header {header}: {e}")
                continue
            
            for row in rows:
                if len(row) != len(header):
                    self.logger.warning(f"Row length {len(row)} does not match header length {len(header)}: {row}")
                    continue
                    
                try:
                    smiles = row[smiles_idx]
                    identifier = row[identifier_idx] if identifier_idx < len(row) else None
                    source = row[source_idx] if source_idx is not None and source_idx < len(row) else database
                    
                    if keep_catalogue and identifier:
                        hits_info.append((smiles, (source, identifier)))
                    else:
                        hits_info.append((smiles, None))
                        
                except (IndexError, TypeError) as e:
                    self.logger.warning(f"Error processing row {row}: {e}")
                    continue
        
        if len(hits_info) == 0:
            self.logger.warning(
                f"No superstructures found for {query_smiles}, returning original query reactant, {query_smiles}.")
            # add query smiles just in case it's not returned
            hits_info.append((query_smiles, None))
        return hits_info

    @staticmethod
    def get_resp_json(url: str,
                      api_key: str = None,
                      data: Dict = None,
                      retries: int = 50,
                      backoff_factor: float = 0.5) -> Optional[Dict] | None:
        """
        Directly get the response json from a request, with retry mechanism for handling 429 status code.
        """
        logger = logging.getLogger(__name__)
        for attempt in range(retries):
            try:
                headers = {'Content-Type': 'application/json'}
                if api_key:
                    headers['X-API-KEY'] = api_key

                response = requests.post(
                    url,
                    headers=headers,
                    data=json.dumps(data),
                )

                if response.status_code in [429, 504]:
                    if attempt < retries - 1:
                        # Calculate wait time using jittered exponential backoff strategy with at most 3 minutes
                        wait_time = backoff_factor * (2 ** attempt)
                        if wait_time > 180:
                            # choose randomly num attempts to wait for
                            wait_time = random.uniform(0, 180)
                        error_type = "Rate limit exceeded" if response.status_code == 429 else "Gateway timeout"
                        logger.warning(f"{error_type}. Waiting for {wait_time} seconds before retrying...")
                        time.sleep(wait_time)
                        continue
                    else:
                        logger.error(
                            f"Max retries exceeded with status code {response.status_code}. Please try again later. "
                            f"{response.status_code}")
                        raise APIRetryLimitExceeded(smiles=data.get('query', None))

                response.raise_for_status()
                return response.json()

            except requests.exceptions.HTTPError as err:
                logger.error(f"HTTP error: {err}")
            except requests.exceptions.ConnectionError as err:
                logger.error(f"Connection error: {err}")
            except requests.exceptions.Timeout as err:
                logger.error(f"Timeout error: {err}")
            except requests.exceptions.RequestException as err:
                logger.error(f"Error: {err}")
                break  # Exit the loop on non-recoverable errors

        return None

    @staticmethod
    def get_search_results(url: str,
                           data: Dict[str, Any],
                           api_key: str = None,
                           page: int = 1,
                           max_pages: int = 900) -> List[Dict] | None:
        """
        Recursively get all pages for the endpoint until reach null next page or
        the max_pages threshold. The default max_pages is set to 900 since the recursion limit is 1000.
        """
        current_limit = sys.getrecursionlimit()
        if max_pages > (current_limit - 100):
            raise ValueError(f"max_pages ({max_pages}) too close to recursion limit ({current_limit})")
        # List where we will gather up all the hits_path.
        all_hits = []
        data = {
            **data,
            'page': page,
        }
        response: Dict | None = Arthor.get_resp_json(url, api_key, data)
        if response is None:  # API error, returning None
            return None
        elif response.get('data') is not None:
            all_hits.extend(response.get('data', []))
        elif response.get('results') is not None:
            all_hits.extend(response.get('results', []))
        # Grab more hits_path if there is a next page supplied.
        if page >= max_pages:
            return all_hits
        next_page = response.get('nextPage', None)
        if next_page is not None:
            next_hits = Arthor.get_search_results(
                url,
                data,
                api_key,
                next_page,
                max_pages
            )
            all_hits.extend(next_hits)
        return all_hits

    def make_query(self):
        """
        This function is used to make a query to the Arthor database.
        """
        pass

    def save_database_search(self):
        """
        This function is used to save the results of the Arthor search as a .csv file.
        """
        pass

    def load_database_search(self):
        """
        This function is used to load the results of the Arthor search from a .csv file.
        """
        pass


__all__ = ['Arthor']
