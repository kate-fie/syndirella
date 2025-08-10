# python

import logging
import numpy as np
import numpy.typing as npt
import pandas as pd
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from scipy.spatial.distance import jaccard
from sklearn.metrics.pairwise import cosine_similarity


def setup_logging(log_level=logging.INFO, log_file=None) -> logging.Logger:
    logger = logging.getLogger()
    logger.setLevel(log_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    return logger


logger = setup_logging()


def calc_cosine_similarity(a: np.array, b: np.array) -> float:
    a_reshaped = a.reshape(1, -1)
    b_reshaped = b.reshape(1, -1)
    return cosine_similarity(a_reshaped, b_reshaped)[0][0]


def calc_jaccard_similarity(a: np.array, b: np.array) -> float:
    return 1 - jaccard(a, b)


def morgan_fp(mol: Chem.rdchem.Mol) -> npt.NDArray[Any]:
    fp1 = AllChem.GetMorganFingerprintAsBitVect(
        mol, useChirality=True, radius=2, nBits=1024
    )
    vec1 = np.array(fp1)
    return vec1


def maccs_fp(mol: Chem.rdchem.Mol) -> npt.NDArray[Any]:
    return np.array(MACCSkeys.GenMACCSKeys(mol))


def get_fp(reaction: str, fp: str = "MACCS", concatenate: bool = True) -> npt.NDArray[Any]:
    reactant_str, product_str = reaction.split(">>")
    reactants = reactant_str.split(".")
    products = product_str.split(".")
    logger.debug(f"Reactants: {reactants}")
    logger.debug(f"Products: {products}")
    reactant_mols = [Chem.MolFromSmarts(reactant) for reactant in reactants]
    product_mols = [Chem.MolFromSmarts(product) for product in products]

    for mol in reactant_mols:
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSSSR(mol)

    for mol in product_mols:
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSSSR(mol)

    if fp.lower() == "maccs":
        reactant_fp = np.sum(np.array([maccs_fp(mol) for mol in reactant_mols]), axis=0)
        product_fp = np.sum(np.array([maccs_fp(mol) for mol in product_mols]), axis=0)
    elif fp.lower() == "morgan":
        reactant_fp = np.sum(
            np.array([morgan_fp(mol) for mol in reactant_mols]), axis=0
        )
        product_fp = np.sum(np.array([morgan_fp(mol) for mol in product_mols]), axis=0)
    else:
        raise KeyError(
            f"Fingerprint {fp} is not yet supported. Choose between MACCS and Morgan"
        )

    if concatenate:
        reaction_fp = np.concatenate((reactant_fp, product_fp))
    else:
        reaction_fp = np.sum((reactant_fp, product_fp), axis=0)

    return reaction_fp


def classify_reaction(reaction_fp: np.array, smirks: pd.DataFrame, similarity_metric: str, top_n: int,
                      threshold: float) -> list[dict]:
    similarities = []

    for _, smirk in smirks.iterrows():
        if similarity_metric == 'cosine':
            sim = calc_cosine_similarity(reaction_fp, smirk['reaction_fp'])
        elif similarity_metric == 'jaccard':
            sim = calc_jaccard_similarity(reaction_fp, smirk['reaction_fp'])
        else:
            raise ValueError(f"Unknown similarity metric: {similarity_metric}")

        similarities.append({
            'name': smirk['name'],
            'similarity': sim
        })

    similarities.sort(key=lambda x: x['similarity'], reverse=True)
    return [s for s in similarities[:top_n] if s['similarity'] >= threshold]
