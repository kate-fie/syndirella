import sys
from typing import Dict, Tuple

from config import config

POSTERA_SEARCHER = None


def postera_similarity_search(mol, thr) -> Dict[str,Tuple[float,Dict]]:
  """

  :param mol:
  :param thr:
  :return: Dictionary of query smiles as keys, and as values a tuple with two elementss. First elem is similarity, second is metadata dict
  """
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
  print(results)
  results = { smi: (similarity, metaData) for similarity, smi, metaData in results if similarity >= thr }
  return results


def postera_superstructure_search(mol, thr) -> Dict[str,Tuple[float,Dict]]:
  if config.POSTERA_MANIFOLD_DIR not in sys.path:
    sys.path.append(config.POSTERA_MANIFOLD_DIR)
  from postera_superstructureSearch import Postera_superstructureSearch
  from rdkit import Chem

  global POSTERA_SEARCHER

  if isinstance(mol, str):
    mol = Chem.MolFromSmiles(mol)

  if POSTERA_SEARCHER is None:
    # TODO: Add this as option to change.
    catalogues = ["enamine_bb"]
    POSTERA_SEARCHER = Postera_superstructureSearch(cache_fname=config.POSTERA_MANIFOLD_SUPERSTRUCTURE_CACHE_FNAME,  verbose=True, catalogues=catalogues,
                                       cache_fname_tags='enamine')
    print('SEARCHING FOR ', Chem.MolToSmiles(mol))
    print('SEARCHING THROUGH ', catalogues)
  results = POSTERA_SEARCHER.search([mol])[0][1]
  print(results)
  # TODO: Change results structure to work without similarity value.
  results = { smi: metaData for smi, metaData in results}
  #results = { smi: (similarity, metaData) for similarity, smi, metaData in results if similarity >= thr }
  return results


