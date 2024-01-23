#!/usr/bin/env python3
"""
syndirella.cobblers_workshop._library.py

This module contains the Library class. This class contains information about the analogue library. It will create the
analogue library from the Reaction object. It will also store the analogue library as a .csv file.
"""

import pandas as pd
from rdkit import Chem
from ._reaction import Reaction
from ._postera import Postera
from typing import (Any, Callable, Union, Iterator, Sequence, List, Dict, Tuple)
import os
from collections import OrderedDict
from rdkit import DataStructs
from rdkit.Chem import AllChem
import glob

class Library:
    """
    This class contains information about the analogue library. It will create the analogue library from the Reaction
    object. It will also store the analogue library as a .csv file.

    """
    def __init__(self, reaction: Reaction, output_dir: str, id: str, num_steps: int, current_step: int):
        self.reaction: Reaction = reaction
        self.output_dir: str = output_dir
        self.id: str = id
        self.num_steps: int = num_steps
        self.current_step: int = current_step
        self.analogues_dataframes: Dict[str: pd.DataFrame] = {}
        self.r1 = None
        self.r2 = None

    def create_library(self):
        """
        This function is used to create the analogue library from the Reaction object.
        """
        for key, value in self.reaction.matched_smarts_to_reactant.items():
            reactant = value[0]
            reactant_smarts = key
            analogue_prefix = value[2]
            df, analogue_columns = self.process_reactant(reactant, reactant_smarts, analogue_prefix)
            df: pd.DataFrame
            analogue_columns: Tuple[str, str]
            self.analogues_dataframes[analogue_prefix] = (df, analogue_columns)

    def process_reactant(self, reactant: Chem.Mol, reactant_smarts: str, analogue_prefix: str) -> pd.DataFrame:
        # search for analogues csv if already created. Perform for all reactants.
        reactant_analogues_path: str = self.check_analogues_csv_exists(analogue_prefix, reactant)
        if reactant_analogues_path is not None:
            analogues: List[str] = self.load_library(reactant_analogues_path, analogue_prefix)
        else: # perform database search
            postera_search = Postera(self.reaction, self.output_dir, search_type="superstructure")
            analogues: List[str] = postera_search.perform_database_search(reactant)
        processed_analogues_df, analogue_columns = (
            self.process_analogues(analogues, reactant_smarts, analogue_prefix))
        processed_analogues_df: pd.DataFrame
        analogue_columns: Tuple[str, str]
        self.save_library(processed_analogues_df, analogue_prefix)
        return processed_analogues_df, analogue_columns

    def process_analogues(self, analogues: List[str], reactant_smarts: str, analogue_prefix: str) -> pd.DataFrame:
        """
        This function puts list of analogues in dataframe and does SMART checking to check if the analogues contains
        the SMARTS pattern of the original reactant and against all other reactants SMARTS.
        """
        analogues: List[str] = self.filter_analogues(analogues, analogue_prefix)
        reactant_smarts_mol: Chem.Mol = Chem.MolFromSmarts(reactant_smarts)
        analogues_mols: List[Chem.Mol] = [Chem.MolFromSmiles(analogue) for analogue in analogues]
        contains_smarts_pattern: List[bool] = \
            [bool(analogue.GetSubstructMatches(reactant_smarts_mol)) for analogue in analogues_mols]
        contains_other_reactant_smarts_pattern, other_reactant_prefix = (
            self.check_analogue_contains_other_reactant_smarts_pattern(analogues_mols, reactant_smarts))
        contains_other_reactant_smarts_pattern: List[bool]
        other_reactant_prefix: str
        analogues_df = pd.DataFrame({f"{analogue_prefix}_smiles": analogues,
                                     f"{analogue_prefix}_mol": analogues_mols,
                                     f"{analogue_prefix}_{self.reaction.reaction_name}": contains_smarts_pattern,
                                     f"{other_reactant_prefix}_{self.reaction.reaction_name}": contains_other_reactant_smarts_pattern})
        return analogues_df, (f"{analogue_prefix}_{self.reaction.reaction_name}", f"{other_reactant_prefix}_{self.reaction.reaction_name}")

    def filter_analogues(self, analogues: List[str], analogue_prefix: str) -> List[str]:
        """
        This function is used to filter out analogues that are not valid.
        """
        mols: List[Chem.Mol] = [Chem.MolFromSmiles(mol) for mol in analogues]
        filtered_mols: List[Chem.Mol] = list(OrderedDict((DataStructs.BitVectToText(AllChem.GetMorganFingerprintAsBitVect(
            mol, 2)), mol) for mol in mols).values())
        filtered_analogues: List[str] = [Chem.MolToSmiles(mol) for mol in filtered_mols]
        num_filtered = len(analogues) - len(filtered_analogues)
        percent_diff = round((num_filtered / len(analogues)) * 100, 2)
        print(f'Removed {num_filtered} replicates ({percent_diff}%) of {analogue_prefix} analogues.')
        return filtered_analogues

    def check_analogue_contains_other_reactant_smarts_pattern(self, analogues_mols: List[Chem.Mol],
                                                              reactant_smarts: str) -> List[bool]:
        """
        This function is used to check if the analogues contains the SMARTS patterns of the other reactant.
        """
        # get other reactant SMARTS pattern
        other_reactant_smarts_pattern = [smarts for smarts in self.reaction.matched_smarts_to_reactant.keys() if not smarts == reactant_smarts][0]
        assert other_reactant_smarts_pattern is not None, "Other reactant SMARTS pattern not found."
        other_reactant_prefix = self.reaction.matched_smarts_to_reactant[other_reactant_smarts_pattern][2]
        other_reactant_smarts_mol: Chem.Mol = Chem.MolFromSmarts(other_reactant_smarts_pattern)
        return [bool(analogue.GetSubstructMatches(other_reactant_smarts_mol)) for analogue in analogues_mols], other_reactant_prefix

    def check_analogue_contains_all_smarts_patterns(self, analogues_mols: List[Chem.Mol]) -> List[bool]:
        """
        This function is used to check if the analogues contains all the SMARTS patterns of the other reactants.
        """
        return NotImplementedError()

    def check_analogues_csv_exists(self, analogue_prefix: str, reactant: Chem.Mol) -> str:
        """
        This function is used to check if the analogue library has already been created and saved as a .csv file.
        """
        if self.current_step == 2 and self.num_steps == 2:
            csv_path: str = self.find_products_csv_name(reactant)
            return csv_path
        os.makedirs(f"{self.output_dir}/extra/", exist_ok=True)
        csv_name = (f"{self.id}_{self.reaction.reaction_name}_"
                    f"{analogue_prefix}_{self.current_step}of{self.num_steps}.csv")
        if os.path.exists(f"{self.output_dir}/extra/{csv_name}"):
            return f"{self.output_dir}/extra/{csv_name}"
        else:
            return None

    def find_products_csv_name(self, reactant: Chem.Mol) -> str:
        """
        This function is used to find the name of the products .csv file by comparing the reactant to the product with
        bit vector similarity.
        """
        product_csvs: List[str] = glob.glob(f"{self.output_dir}/extra/"
                                            f"{self.id}_*products*.csv")
        for i, path in enumerate(product_csvs):
            df = pd.read_csv(path)
            product = df["smiles"].tolist()[0]
            product_mol = Chem.MolFromSmiles(product)
            similarity_score = DataStructs.FingerprintSimilarity(
                AllChem.GetMorganFingerprintAsBitVect(reactant, 2),
                AllChem.GetMorganFingerprintAsBitVect(product_mol, 2))
            if similarity_score == 1:
                return path

    def save_library(self, df: pd.DataFrame, analogue_prefix: str):
        """
        This function is used to save the analogue library as a .csv file in output_dir/extra/.
        """
        os.makedirs(f"{self.output_dir}/extra/", exist_ok=True)
        csv_name = f"{self.id}_{self.reaction.reaction_name}_{analogue_prefix}_{self.current_step}of{self.num_steps}.csv"
        df.to_csv(f"{self.output_dir}/extra/{csv_name}")

    def load_library(self, reactant_analogues_path: str, analogue_prefix: str) -> List[str]:
        """
        This function is used to load the analogue library from a .csv file.
        """
        df = pd.read_csv(reactant_analogues_path)
        # find column with analogue smiles
        try:
            analogues = df["smiles"].tolist()
            return analogues
        except KeyError:
            print(f"No 'smiles' column found in {reactant_analogues_path}. "
                  f"Looking for {analogue_prefix}_smiles column.")
        try:
            analogues = df[f"{analogue_prefix}_smiles"].tolist()
            return analogues
        except KeyError:
            print(f"No {analogue_prefix}_smiles column found in {reactant_analogues_path}. Stopping...")