#!venv/bin/env python3
"""
slipper_fitter/CobblersWorkshop.py

This module contains the SlipperFitter class.
"""
from typing import (List, Dict, Tuple, Union, Optional)
from fragmenstein import Laboratory, Wictor, Igor
import pandas as pd
from fragmenstein.laboratory.validator import place_input_validator
from rdkit import Chem
from rdkit.Chem import rdinchi
import os, logging
import time
from syndirella.slipper import intra_geometry, flatness
import datetime

class SlipperFitter:
    """
    This class is instantiated to place all the products in the template using hits_path
    """

    def __init__(self,
                 template_path: str,
                 hits_path: str,
                 hits_names: List[str],
                 output_dir: str):
        self.template_path: str = template_path
        self.hits_path: str = hits_path
        self.hits_names: List[str] = hits_names
        self.output_dir: str = output_dir
        self.num_atom_diff_limit: int = 10
        self.final_products_csv_path: str = None
        self.final_products_pkl_path: str = None
        self.batch_num: int = None
        self.final_products: pd.DataFrame = None
        self.placements: pd.DataFrame = None
        self.merged_placements: pd.DataFrame = None
        self.output_path = None
        # Placements variables set
        self.n_cores: int = 8
        self.timeout: int = 240

    def fit_products(self):
        """
        Main entry to the SlipperFitter class. This function is used to fit the products to the final library.
        """
        print('Fitting products, warning: it is expected that the base compound has already passed minimisation and '
              'intramolecular checks.')
        input_df: pd.DataFrame = self.prep_products()
        placements = self.place_products(input_df)
        self.placements = self.check_intra_geometry(placements)
        self.edit_placements()
        self.merged_placements = self.merge_placements()
        return self.merged_placements

    def check_base(self, base: Chem.Mol) -> bool:
        """
        Checks if the base can be minimised (no specific stereoisomer) and passes intermolecular checks.
        If it cannot be minimised after 3 attempts, returns False.
        """
        input_df: pd.DataFrame = self._prep_base_input_df(base)
        id = generate_inchi_ID(Chem.MolToSmiles(base))
        output_path: str = os.path.join(self.output_dir, f'{id}-base-check')
        lab: Laboratory = self.setup_Fragmenstein(output_path)
        num_attempts = 2
        for attempt in range(1, num_attempts+1):
            base_placed: Chem.Mol = self._place_base(lab, input_df)  # None if not minimised
            if base_placed is not None:
                geometries: Dict = intra_geometry.check_geometry(base_placed,
                                                                 threshold_clash=0.4) # increasing threshold for internal clash
                flat_results: Dict = flatness.check_flatness(base_placed)
                if self._check_intra_geom_flatness_results(geometries=geometries, flat_results=flat_results):
                    return True
                    if len(os.listdir(f'{output_path}/base-check')) > 0: # last resort just check if there are files
                            print('Base could be minimised and passed intramolecular checks!')
                            return True
                    else:
                        print(f'Base could not be minimised. Attempt {attempt} of {num_attempts}.')
                else:
                    print(f'Base could not pass intramolecular checks. Attempt {attempt} of {num_attempts}.')
            else:
                print(f'Base could not be minimised. Attempt {attempt} of {num_attempts}.')
        return False

    def _check_intra_geom_flatness_results(self,
                                           geometries: Dict,
                                           flat_results: Dict) -> bool:
        """
        This function checks the intramolecular geometry and flatness results and returns True if all are True.
        """
        if not geometries['results']['bond_lengths_within_bounds']:
            print('Mol could not pass bond length checks.')
            return False
        if not geometries['results']['bond_angles_within_bounds']:
            print('Mol could not pass bond angle checks.')
            return False
        if not geometries['results']['no_internal_clash']:
            print('Mol could not pass internal clash checks.')
            return False
        if not flat_results['results']['flatness_passes']:
            print('Mol could not pass flatness checks.')
            return False
        return True

    def _prep_base_input_df(self, base: Chem.Mol) -> pd.DataFrame:
        """
        Prepares input dataframe to Fragmenstein for the base compound.
        """
        base_df: pd.DataFrame = pd.DataFrame({'name': ['base_check'],
                                              'smiles': [Chem.MolToSmiles(base)]})
        base_df = self.add_hits(base_df)
        base_df = base_df.astype({'name': str, 'smiles': str})
        return base_df

    def _get_base(self, input_df: pd.DataFrame) -> Chem.Mol:
        """
        Get base compound as mol object from the input_df.
        """
        # get flat base smiles
        base_smiles: str = input_df.loc['base' in input_df['name'], 'smiles'].values[0]
        base_mol: Chem.Mol = Chem.MolFromSmiles(base_smiles)
        return base_mol

    def prep_products(self) -> pd.DataFrame:
        """
        This function is used to prepare the inputs for Fragmenstein.
        """
        # input_df is a copy of final_products but removing duplicates of names and other metadata from synthesizer step
        input_df: pd.DataFrame = self._prepare_input_df()
        input_df = self.add_hits(input_df)
        # add typing to columns
        input_df = input_df.astype({'name': str, 'smiles': str})
        # Check if there are any duplicates in the input_df
        if input_df.duplicated(subset='name').any():
            raise ValueError('There are duplicates of names in the product dataframe to place.')
        return input_df

    def _prepare_input_df(self) -> pd.DataFrame:
        """
        Creates the input dataframe for the placements. Will remove duplicates of names and other metadata from
        synthesizer step.
        """
        input_df: pd.DataFrame = self.final_products.copy(deep=True)
        # drop duplicates of names
        input_df = input_df.drop_duplicates(subset='name')
        # cut rows with number of atoms difference greater than num_atom_diff_limit
        print()
        print(f'Cutting products with number of atoms difference greater than {self.num_atom_diff_limit}.')
        input_df = input_df[input_df['num_atom_diff'] <= self.num_atom_diff_limit]
        self._print_diff(self.final_products, input_df, verb='Kept')
        # drop columns that are not needed
        input_df = input_df[['name', 'smiles']]
        # place number in batch
        if len(input_df) > self.batch_num:
            print(f'Even after cutting products, the number of products is over {self.batch_num}.'
                  f' Only placing the top {self.batch_num} products.')
            input_df = input_df.iloc[:self.batch_num]
        self._print_diff(self.final_products, input_df, verb='Placing')
        return input_df

    def _print_diff(self,
                    orig_df: pd.DataFrame,
                    input_df: pd.DataFrame,
                    verb: str = None):
        """
        This function is used to print the difference between the original number of analogues and the number of
        valid analogues.
        """
        assert len(input_df) <= len(orig_df), ("Problem with finding unique analogues. There are more than were in "
                                               "the original list of analogues")
        num_placed = len(input_df)
        percent = round(((len(input_df) / len(orig_df)) * 100), 2)
        print(f'{verb} {len(input_df)} ({percent}%) unique analogues out of {len(orig_df)} analogues.')

    def add_hits(self, input_df: pd.DataFrame) -> pd.DataFrame:
        """
        This function adds the hits_path as mol objects to input_df['hits_path'].
        """
        # load hits_path either from mol or sdf
        if os.path.splitext(self.hits_path)[1] == '.mol':
            print('This is a mol file')
            hits: List[Chem.Mol] = [Chem.MolFromMolFile(self.hits_path.strip())]
        else:
            with Chem.SDMolSupplier(self.hits_path.strip()) as sd:
                hits: List[Chem.Mol] = list(sd)
        # Find which hits are in the hit_names
        hits = [
            hit for hit in hits
            if any(hit_name in hit.GetProp('_Name') for hit_name in self.hits_names)
        ]
        # add hits column with contents of hits_names as space separated string
        input_df['hits'] = input_df.apply(lambda row: hits, axis=1)
        return input_df

    def _place_base(self,
                    lab: Laboratory,
                    input_df: pd.DataFrame) -> Chem.Mol:
        """
        Places the base compound, returns the mol object of the base compound if successful else None.
        """
        placements: pd.DataFrame = lab.place(place_input_validator(input_df), n_cores=self.n_cores,
                                             timeout=self.timeout)
        try:
            base: Chem.Mol = placements.at[0, 'minimized_mol']
        except Exception:
            base = None
        return base

    def place_products(self, input_df: pd.DataFrame) -> pd.DataFrame:
        """
        This function places products using Fragmenstein.
        """
        start_time = time.time()  # Start timing
        # set up Wictor
        self.output_path: str = self.output_dir + '/output'
        lab: Laboratory = self.setup_Fragmenstein(self.output_path)
        placements: pd.DataFrame = lab.place(place_input_validator(input_df),
                                             n_cores=self.n_cores,
                                             timeout=self.timeout)
        end_time = time.time()  # End timing
        elapsed_time = end_time - start_time  # Calculate elapsed time
        print(f"Placing {len(input_df)} run time: {datetime.timedelta(seconds=elapsed_time)}")
        return placements

    def check_intra_geometry(self,
                             placements: pd.DataFrame) -> pd.DataFrame:
        """
        Checks the intramolecular geometry and flatness of double bonds and aromatic rings
        of each mol object using PoseBuster's implemented checks.

        Checks:
        Bond lengths
            The bond lengths in the input molecule are within 0.75 of the lower and 1.25 of the upper bounds determined by distance geometry
        Bond angles
        	The angles in the input molecule are within 0.75 of the lower and 1.25 of the upper bounds determined by distance geometry
        Planar aromatic rings
            All atoms in aromatic rings with 5 or 6 members are within 0.25 Å of the closest shared plane
        Planar double bonds
            The two carbons of aliphatic carbon–carbon double bonds and their four neighbours are within 0.25 Å of the closest shared plane
        Internal steric clash
            The interatomic distance between pairs of non-covalently bound atoms is above 0.7 of the lower bound determined by distance geometry
        """
        # get list of mol objects
        mols: List[Chem.Mol] = placements['minimized_mol'].tolist()
        assert len(mols) == len(placements), "There are missing mol objects in the placements."
        intra_geometry_results: List[bool] = []
        flatness_results: List[bool] = []
        for mol in mols:
            if mol is None:
                intra_geometry_results.append(None)
                flatness_results.append(None)
            else:
                geometries: Dict = intra_geometry.check_geometry(mol, threshold_clash=0.4)  # increasing threshold for internal clash
                flat_results: Dict = flatness.check_flatness(mol)
                intra_geometry_results.append(self._check_intra_geom_flatness_results(geometries=geometries, flat_results=flat_results))
        # get list of pass output from intra_geometry for each mol
        placements['intra_geometry_pass'] = intra_geometry_results
        return placements

    def setup_Fragmenstein(self,
                           output_path: str) -> Laboratory:
        """
        This function sets up Fragmenstein to run.
        """
        # Using Wictor to place just by RDKit minimisation
        # make output directory
        os.makedirs(output_path, exist_ok=True)
        Wictor.work_path = output_path
        os.chdir(output_path)  # this does the trick
        Wictor.monster_throw_on_discard = True  # stop this merger if a fragment cannot be used.
        Wictor.monster_joining_cutoff = 5  # Å
        Wictor.quick_reanimation = False  # for the impatient
        Wictor.error_to_catch = Exception  # stop the whole laboratory otherwise
        Wictor.enable_stdout(logging.CRITICAL)
        # Wictor.enable_logfile(os.path.join(self.output_path, f'fragmenstein.log'), logging.ERROR)
        Laboratory.Victor = Wictor
        Igor.init_pyrosetta()
        with open(self.template_path) as fh:
            pdbblock: str = fh.read()
        lab = Laboratory(pdbblock=pdbblock,
                         covalent_resi=None,
                         run_plip=False)
        return lab

    def edit_placements(self):
        """
        This function edits the placements.
        """
        self.placements.loc[(self.placements['∆∆G'] < 0) & (self.placements['comRMSD'] < 2), 'outcome'] = 'acceptable'
        self.placements.loc[(self.placements['∆∆G'] > -1) &
                            (self.placements.outcome == 'acceptable'), 'outcome'] = 'weak'
        self.placements['unminimized_mol'] = self.placements.unminimized_mol.fillna(Chem.Mol())
        self.fix_intxns()
        if self.placements.outcome.value_counts().get('acceptable') is None:
            num_success = 0
            percent_success = 0
        else:
            num_success = self.placements.outcome.value_counts()['acceptable']
            percent_success = round((self.placements.outcome.value_counts()['acceptable'] / len(self.placements) * 100)
                                    , 2)
        print(f'{num_success} ({percent_success}%) successful placements '
              f'where ∆∆G < -1 and RMSD < 2 Å.')
        print()

    def fix_intxns(self):
        intxn_names: List[tuple] = [c for c in self.placements.columns if isinstance(c, tuple)]
        for intxn_name in intxn_names:
            self.placements[intxn_name] = self.placements[intxn_name].fillna(0).astype(int)

    def merge_placements(self) -> pd.DataFrame:
        """
        This function merges the fragmenstein output with products csv.
        """
        # Define the columns to drop
        columns_to_drop = ['smiles', 'binary_hits', 'mode', 'runtime', 'disregarded', 'unmin_binary', 'min_binary',
                           'hit_binaries', 'unminimized_mol', 'minimized_mol', 'hit_mols', 'hit_names']
        # Filter out the columns that do not exist in merged_placements
        columns_to_drop = [col for col in columns_to_drop if col in self.placements.columns]
        # Drop the columns if they exist
        self.placements = self.placements.drop(columns=columns_to_drop)
        # merge on name
        merged_placements = pd.merge(self.final_products, self.placements, on='name', how='left')
        return merged_placements

    def save_placements(self):
        """
        This function saves the placements as a .pkl.gz and .csv file.
        """
        self.placements.to_pickle(os.path.join(self.output_dir, 'fragmenstein_placements.pkl.gz'))
        self.placements.to_csv(os.path.join(self.output_dir, 'fragmenstein_placements.csv'))
        # get final name same as final products csv
        merged_placements_csv_path = self.final_products_csv_path.split('.')[0] + '_placements.csv'
        merged_placements_pkl_path = self.final_products_pkl_path.split('.')[0] + '_placements.pkl.gz'
        self.merged_placements.to_csv(merged_placements_csv_path)
        self.merged_placements.to_pickle(merged_placements_pkl_path)


def generate_inchi_ID(smiles: str) -> str:
    """
    This function is used to generate a unique id for the route just using the product.
    """
    assert Chem.MolFromSmiles(smiles), f"Could not create a molecule from the smiles {smiles}."
    ID = rdinchi.MolToInchi(Chem.MolFromSmiles(smiles))
    id = rdinchi.InchiToInchiKey(ID[0])
    return id
