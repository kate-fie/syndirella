#!venv/bin/env python3
"""
syndirella.slipper.Slipper.py

This module contains the Slipper class. A slipper in this metaphor is the set of molecules that is the
product of a reaction.
"""
from typing import List, Dict, Optional
import pandas as pd
from syndirella.slipper.slipper_synthesizer.SlipperSynthesizer import SlipperSynthesizer
from syndirella.cobblers_workshop.Library import Library
from syndirella.slipper.SlipperFitter import SlipperFitter
from syndirella.slipper._placement_data import get_placement_data
import os, shutil
import glob
from rdkit import Chem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from rdkit.Chem import AllChem, inchi


class Slipper:
    """
    This class is instantiated to represent all products for a step in a route.
    """
    def __init__(self,
                 *,
                 library: Library,
                 template: str = None,
                 hits: str = None,
                 hits_names: List[str] = None,
                 batch_num: int = None,
                 atoms_ids_expansion: dict = None,
                 additional_info: dict = None):
        self.products: pd.DataFrame = None
        self.library: Library = library
        self.output_dir: str = library.output_dir
        self.final_products_pkl_path: str = None
        self.final_products_csv_path: str = None
        # need Fragmenstein information
        self.template: str = template # path to pdb file
        self.hits: str = hits # path to .sdf or .mol file
        self.hits_names: List[str] = hits_names # name of fragments
        self.batch_num: int = batch_num
        self.atoms_ids_expansion: dict = atoms_ids_expansion
        self.placements: pd.DataFrame = None
        self.output_path: str = None
        self.additional_info: dict = additional_info

    def get_products(self) -> pd.DataFrame and str:
        """
        Main entry to the Slipper class. This function is used to get the products from the final library.
        """
        slipper_synth = SlipperSynthesizer(self.library,
                                           self.output_dir,
                                           self.atoms_ids_expansion,
                                           self.additional_info)
        self.products: pd.DataFrame = slipper_synth.get_products()
        uuid = slipper_synth.uuid
        if self.atoms_ids_expansion is not None:
            slipper_synth.label_products()
        slipper_synth.save_products()
        self.final_products_pkl_path: str = slipper_synth.final_products_pkl_path
        self.final_products_csv_path: str = slipper_synth.final_products_csv_path
        return self.products, uuid

    def place_products(self):
        """
        This function is used to place the products with Fragmenstein.
        """
        slipper_fitter = SlipperFitter(self.template,
                                       self.hits,
                                       self.hits_names,
                                       self.output_dir)
        slipper_fitter.final_products = self.products
        slipper_fitter.batch_num = self.batch_num
        slipper_fitter.final_products_pkl_path = self.final_products_pkl_path
        slipper_fitter.final_products_csv_path = self.final_products_csv_path
        self.placements: pd.DataFrame = slipper_fitter.fit_products()
        slipper_fitter.save_placements()
        self.output_path: str = slipper_fitter.output_path
        return self.placements

    def write_products_to_hippo(self,
                                uuid) -> str:
        """
        Writes a dataframe that contains the values needed for HIPPO db input.

        Open questions:
        - I think it should be written at the end of each reaction step.
        - Should be checked if the dataframe already exists and you're adding to it

        Returns:
            path: str : the path to the saved dataframe
        """
        assert self.placements is not None, "Placements need to be run first before writing HIPPO output."
        # cut placements to those that were placed by batch_num
        placements: pd.DataFrame = self.placements.iloc[:self.batch_num]

        hippo_path: str = self.output_dir + f'/{self.library.id}_{uuid}_to_hippo.pkl.gz'
        # get all products dfs in /extra
        products_files: List[str] = glob.glob(f"{self.output_dir}/extra/*products*.pkl*")
        product_dfs: Dict[int, pd.DataFrame] = self._load_products_dfs(products_files)
        # make HIPPO output dataframe of these specific products
        hippo_df = self._structure_products_for_hippo(placements_df=placements,
                                                      product_dfs=product_dfs)
        # load file if it already exists
        if os.path.isfile(hippo_path):
            previous_hippo = pd.read_pickle(hippo_path)
            hippo_df = pd.concat([previous_hippo, hippo_df], axis=1, ignore_index=True)
        # save the dataframe
        hippo_df.to_pickle(hippo_path)
        print(f"Saved HIPPO output to {hippo_path}")
        return hippo_path

    def _load_products_dfs(self, products_files: List[str]) -> Dict[int, pd.DataFrame]:
        """
        Load the products dataframes from the files in the /extra directory and putting into dict where key is step.
        """
        product_dfs: Dict[int, pd.DataFrame] = {}
        for file in products_files:
            df = pd.read_pickle(file)
            step = df['step'].iloc[0]
            product_dfs[step] = df
        assert len(product_dfs) == self.library.num_steps-1, "Not all steps have products dataframes."
        return product_dfs

    def _structure_products_for_hippo(self,
                                      placements_df: pd.DataFrame,
                                      product_dfs: Dict[int, pd.DataFrame]) -> pd.DataFrame:
        """
        Structures the placements or products df for HIPPO output.
        """
        hippo_dfs: Dict[int, pd.DataFrame] = {}
        for step, products_df in product_dfs.items():
            hippo_df_step: pd.DataFrame = self._structure_step_for_hippo(step, products_df)
            hippo_dfs[step] = hippo_df_step
        hippo_dfs[self.library.num_steps] = self._structure_step_for_hippo(self.library.num_steps, placements_df)
        hippo_df = self._put_hippo_dfs_together(hippo_dfs)
        return hippo_df

    def calculate_inchi_similarity(self,
                                   smiles1: str,
                                   smiles2: str) -> int:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        # Check if either molecule failed to be created from the SMILES
        if mol1 is None or mol2 is None:
            return 0  # Indicating an error or non-comparability
        inchi1 = inchi.MolToInchi(mol1)
        inchi2 = inchi.MolToInchi(mol2)
        # Directly compare the InChI strings for equality
        similarity = 1 if inchi1 == inchi2 else 0
        return similarity

    def calculate_tanimoto(self,
                           smiles1,
                           smiles2):
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        if mol1 is None or mol2 is None:
            return 0
        fp1 = AllChem.GetMorganFingerprint(mol1, 2)
        fp2 = AllChem.GetMorganFingerprint(mol2, 2)
        return TanimotoSimilarity(fp1, fp2)

    # Function to find matches
    def find_matches(self,
                     row, # row of current step
                     step,
                     df_previous_step: pd.DataFrame) -> List[str]:
        product_matches = None  # Store matching product names
        for _, product_row in df_previous_step.iterrows():
            similarity = self.calculate_inchi_similarity(row[f'{step + 1}_r{row[f"{step + 1}_r_previous_product"]}_smiles'],
                                                         product_row[f'{step}_product_smiles'])
            if similarity == 1:
                product_matches = product_row[f'{step}_product_name']
                break
        return product_matches


    def _put_hippo_dfs_together(self,
                                hippo_dfs: Dict[int, pd.DataFrame]) -> pd.DataFrame:
        """
        Puts the HIPPO dataframes together by matching on each reaction product to the correct previous step's reactant.
        """
        # get the last step's dataframe
        hippo_df_step_last = hippo_dfs[self.library.num_steps]
        for step in range(1, self.library.num_steps)[::-1]: # iterate through the steps in reverse
            # get the step's dataframe
            hippo_df_stepx = hippo_dfs[step]
            # Find the matching product names for the reactant in this new step
            hippo_df_step_last[f'{step}_product_name'] = hippo_df_step_last.apply(self.find_matches,
                                                                                      step=step,
                                                                                      df_previous_step=hippo_df_stepx,
                                                                                      axis=1)
            # What happens if there are null product names?... Still keep row with null product name
            # Join on the product names, has to be right merge because we only care about the products from the last step
            result_df = pd.merge(hippo_df_stepx,
                                 hippo_df_step_last,
                                 left_on=f'{step}_product_name',
                                 right_on=f'{step}_product_name',
                                 how='right')
            # Update the last step dataframe
            hippo_df_step_last = result_df
        # add base compound smiles as first column, get from Inchi Key
        base_compound_smiles: str = Chem.MolToSmiles(self.library.reaction.product)
        hippo_df_step_last.insert(0, 'base_compound_smiles', base_compound_smiles)
        return hippo_df_step_last


    def _structure_step_for_hippo(self,
                                  step: int,
                                  products_df: pd.DataFrame) -> pd.DataFrame:
        """
        Structures the products df for HIPPO output.
        """
        # cut down the dataframe to those with num_atom_diff <= 15
        products_df = products_df[products_df['num_atom_diff'] <= 15]
        reaction: str = products_df['reaction'].iloc[0]
        r1_smiles: List[str] = products_df['r1_smiles'].tolist()
        r1_is_previous_product: bool = products_df['r1_is_previous_product'].iloc[0]
        try:
            r2_smiles: List[str] = products_df['r2_smiles'].tolist()
            r2_is_previous_product: bool = products_df['r2_is_previous_product'].iloc[0]
        except KeyError:
            r2_smiles: List[str] = [''] * len(products_df)
            r2_is_previous_product: bool = None
        # find which reactants were products of the previous step
        r_previous_product: int = self.which_reactant_was_previous_product(r1_is_previous_product,
                                                                           r2_is_previous_product)
        product_smiles: List[str] = products_df['smiles'].tolist()
        product_names: List[str] = products_df['name'].tolist()
        try:
            flags: List[str] = products_df['flag'].tolist()
            # replace nan with None
            flags = [None if pd.isna(flag) else flag for flag in flags]
        except KeyError:
            flags: List[str] = [None] * len(products_df)
        # variables for the last step
        if step == self.library.num_steps:
            num_atom_diff: List[int] = products_df['num_atom_diff'].tolist()
            stereoisomer: List[str] = products_df['stereoisomer'].tolist()
            error: List[str] = products_df['error'].tolist()
            delta_delta_G: List[float] = products_df['∆∆G'].tolist()
            delta_G_bound: List[float] = products_df['∆G_bound'].tolist()
            delta_G_unbound: List[float] = products_df['∆G_unbound'].tolist()
            comRMSD: List[float] = products_df['comRMSD'].tolist()
            regarded: List[bool] = products_df['regarded'].tolist()
            intra_geometry_pass: List[bool] = products_df['intra_geometry_pass'].tolist()
            path_to_mol: List[Optional[str]] = [
                f'{self.output_dir}/output/{name}/{name}.minimised.mol' if not pd.isna(ddg) else None
                for name, ddg in zip(product_names, products_df['∆∆G'])
            ]
        else:
            num_atom_diff: List[None] = [None] * len(products_df)
        # make HIPPO output dataframe
        hippo_df_step = pd.DataFrame({f'{step}_reaction': reaction,
                                     f'{step}_r1_smiles': r1_smiles,
                                     f'{step}_r2_smiles': r2_smiles,
                                     f'{step}_r_previous_product': r_previous_product,
                                     f'{step}_product_smiles': product_smiles,
                                     f'{step}_product_name': product_names,
                                     f'{step}_num_atom_diff': num_atom_diff,
                                     f'{step}_flag': flags})
        if step == self.library.num_steps:
            hippo_df_step[f'{step}_stereoisomer'] = stereoisomer
            hippo_df_step[f'error'] = error
            hippo_df_step[f'∆∆G'] = delta_delta_G
            hippo_df_step[f'∆G_bound'] = delta_G_bound
            hippo_df_step[f'∆G_unbound'] = delta_G_unbound
            hippo_df_step[f'comRMSD'] = comRMSD
            hippo_df_step[f'regarded'] = regarded
            hippo_df_step[f'path_to_mol'] = path_to_mol
            hippo_df_step[f'intra_geometry_pass'] = intra_geometry_pass
        return hippo_df_step

    def which_reactant_was_previous_product(self,
                                            r1_is_previous_product: bool,
                                            r2_is_previous_product: bool) -> int:
        """
        Determine which reactant was the product of the previous step.
        """
        if r1_is_previous_product:
            assert not r2_is_previous_product, "Both reactants cannot be products of the previous step."
            return 1
        elif r2_is_previous_product:
            return 2
        else:
            return None

    def _delete_file_or_directory(self, path):
        """
        Delete a file or directory at the given path.
        """
        try:
            if os.path.isfile(path) or os.path.islink(path):
                os.unlink(path)
                #print(f"Deleted file: {path}")
            elif os.path.isdir(path):
                shutil.rmtree(path)
                #print(f"Deleted directory: {path}")
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (path, e))

    def _should_delete_file(self, file, suffixes_to_keep):
        """
        Determine if a file should be deleted based on its suffix.
        """
        return not any(file.endswith(suffix) for suffix in suffixes_to_keep)

    def clean_up_placements(self):
        """
        This function is used to remove extra files that are generated by Fragmenstein.
        """
        print("Cleaning up placement directory...")
        suffixes_to_keep = ['.minimised.json', '.minimised.mol', '.csv', '.pkl', '.pkl.gz']
        for root, dirs, files in os.walk(self.output_path):
            for file in files:
                if self._should_delete_file(file, suffixes_to_keep):
                    file_path = os.path.join(root, file)
                    self._delete_file_or_directory(file_path)

    def get_placements_df(self) -> pd.DataFrame:
        """
        This function is used to get the placements dataframe which is the merged df with the products.
        """
        if self.placements is not None:
            return self.placements
        # otherwise you're gonna have to build it from scratch by going through each file in the output dir
        # get products df
        assert self.products is not None, "You need to set the products df first."
        self.placements = get_placement_data(self.products, self.output_path, self.output_dir)
        return self.placements


