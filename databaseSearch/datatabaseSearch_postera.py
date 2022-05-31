import sys
from typing import Dict, Tuple

from config import config

POSTERA_SEARCHER = None


def postera_similarity_search(mol, thr) -> Dict[str,Tuple[float,Dict]]:
  if config.POSTERA_MANIFOLD_DIR not in sys.path:
    sys.path.append(config.POSTERA_MANIFOLD_DIR)
  from postera_similaritySearch import Postera_similaritySearch
  from rdkit import Chem

  global POSTERA_SEARCHER

  if isinstance(mol, str):
    mol = Chem.MolFromSmiles(mol)

  if POSTERA_SEARCHER is None:
    POSTERA_SEARCHER = Postera_similaritySearch(cache_fname=config.POSTERA_MANIFOLD_CACHE_FNAME,  verbose=False)
  results = POSTERA_SEARCHER.search([mol])[0][1]
  results = { smi: (similarity, metaData) for similarity, smi, metaData in results if similarity >= thr }
  return results


def postera_substructure_search(mol, thr):
  raise NotImplementedError()


