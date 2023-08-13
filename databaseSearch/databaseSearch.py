#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 15:10:15 2021

@author: ruben
"""
import importlib
from collections import OrderedDict
from typing import Tuple, Dict

import sys
sys.path.append('/Users/kate_fieseler/PycharmProjects/chemUtils')
sys.path.append('/Users/kate_fieseler/PycharmProjects/retrievesynthesizable/databaseSearch')
import chemUtils
from chemUtils.fragmentation import takeLargestSubFragement
from rdkit import Chem

from config import config


def clean_results(raw_searchResults: Dict[str,Tuple[float,Dict]]) -> Dict[str,Tuple[float,Dict]]:
  results = OrderedDict(zip(map(Chem.MolToSmiles, map(takeLargestSubFragement,
                                                      map(Chem.MolFromSmiles, raw_searchResults)),
                                ),
                            raw_searchResults.values()))
  return results


p, m = config.DATABASE_SEARCH_FUN.rsplit('.', 1)

mod = importlib.import_module(p)
search_function = getattr(mod, m)
search_function_substructure = getattr(mod, m.replace("similarity", "superstructure"))

def perform_database_search(mol: Chem.Mol, substructure:bool, thr:float) -> Dict[str,Tuple[float,Dict]]:
  """

  :param mol: query molecule for database search
  :param thr: search metric to rule out compounds. Compounds are ruled out if mol.val < thr. Generally refers to fingerprint similarity
  :return: A dictionary with smiles as keys, and as values a tuple with two elementss. First elem is similarity, second is metadata dict
  """
  if substructure:
      raw_results = search_function_substructure(mol, thr)
  else:
      raw_results = search_function(mol, thr)
  return clean_results(raw_results)


