#!/usr/bin/env python3
"""
Hippo.py

This module contains the functionality for an Hippo search.
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
# from hippo import HIPPO  # HIPPO dependency removed


class Hippo(DatabaseSearch):
    """
    This class contains information about the Hippo search. It will perform the Hippo search using the
    perform_database_search function. It will also store the results of the Hippo search as a .csv file.
    """

    def __init__(self, db_path: str):
        super().__init__()

        self.logger = logging.getLogger(f"{__name__}")

        self.db_path = db_path        
        # self.animal = HIPPO("reference_db", db_path)  # HIPPO dependency removed
        self.animal = None  # Placeholder for HIPPO functionality

    def perform_database_search(self,
                                reactant: Chem.Mol,
                                reaction_name: str,
                                search_type: str = "superstructure") -> List[str] | None:
        """
        This function is used to perform the Hippo search using the database_search_function.
        """
        if search_type != "superstructure":
            raise NotImplementedError(f"Search type '{search_type}' is not implemented. Only 'superstructure' search is supported.")

        # 1. Get additional similar reactant if reaction is one with additional reactants
        reactant_filters = fairy.load_reactant_filters()
        reactants: List[str] = fairy.find_similar_reactants(reactant=reactant,
                                                            reaction_name=reaction_name,
                                                            reactant_filters=reactant_filters)
        # 2. Perform the search for all
        hits_all: List[Tuple[str, Tuple[str, int]]] = []
        for smiles in reactants:
            if search_type == "superstructure":
                hits: List[Tuple[str, Tuple[str, int]]] = self.perform_superstructure_search(smiles)
                self.logger.info(f'Found {len(hits)} hits for {smiles} before filtering.')
                hits_all.extend(hits)
        filtered_hits: List[str] = fairy.filter_molecules(hits=hits_all)
        return filtered_hits

    def perform_superstructure_search(self,
                                      smiles: str) -> List[Tuple[str, Tuple[str, int]]]:
        """
        This function is used to perform the Hippo superstructure search.
        """

        self.logger.info(f"Running superstructure search for {smiles}.")
        if not isinstance(smiles, str):
            self.logger.error("Smiles must be a string.")
            raise TypeError("Smiles must be a string.")

        if self.animal is None:
            self.logger.warning("HIPPO functionality is not available. Returning empty results.")
            hits = []
        else:
            hits = self.animal.db.query_substructure(smiles, none="quiet")
        self.logger.info(f"#results = {len(hits or [])}")

        hits_info: List[Tuple[str, Tuple[str, int]]] = self.structure_output(hits, query_smiles=smiles)

        return hits_info


    def structure_output(self, 
                         hits: List[Dict] | None, 
                         query_smiles: str) -> List[Tuple[str, Tuple[str, int]]]:
        """
        Formats output into a list of tuples with smiles and catalogue.
        """

        hits_info: List[Tuple[str, Tuple[str, int]]] = []
        if hits is None:
            hits = []

        for compound in hits:
            smiles = compound.smiles
            identifier = compound.id
            hits_info.append((smiles, (self.db_path, identifier)))
                                    
        if len(hits_info) == 0:
            self.logger.warning(
                f"No superstructures found for {query_smiles}, returning original query reactant, {query_smiles}.")
            hits_info.append((query_smiles, (self.db_path, None)))
        return hits_info

    def make_query(self):
        """
        This function is used to make a query to the Hippo database.
        """
        pass

    def save_database_search(self):
        """
        This function is used to save the results of the Hippo search as a .csv file.
        """
        pass

    def load_database_search(self):
        """
        This function is used to load the results of the Hippo search from a .csv file.
        """
        pass


__all__ = ['Hippo']
