import sys

from config import config


def chembl_similarity_search(wholeMol, thr):
  from chembl_webresource_client.new_client import new_client
  from rdkit import Chem

  if isinstance(wholeMol, str):
    smi = wholeMol
  else:
    smi = Chem.MolToSmiles(wholeMol)

  if isinstance(thr, float):
    thr = int(round(thr)*100)

  similarity_query = new_client.similarity  #substructure = new_client.substructure
  results = similarity_query.filter(smiles=smi, similarity=thr).only(['molecule_chembl_id', 'similarity', 'molecule_structures'])
  results = { res["molecule_structures"]["canonical_smiles"]: (res['molecule_chembl_id'], res['similarity'])  for res in results }
  assert len(results) > 0, "Error, database search chembl_similarity_search failed!"
  return results