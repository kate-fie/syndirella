#!venv/bin/env python3
"""
syndirella.slipper._base.py

This module contains the Slipper class. A slipper in this metaphor is the set of molecules that is the
product of a reaction.
"""
from typing import (List, Dict, Tuple, Union, Optional)
import pandas as pd
from syndirella.slipper.slipper_synthesizer._base import SlipperSynthesizer
from syndirella.cobblers_workshop._library import Library
from syndirella.slipper.slipper_fitter._base import SlipperFitter
from syndirella.slipper._placement_data import get_placement_data, get_delta_delta_G, get_bound_unbound, make_success_csv_row
import os, shutil

class Slipper:
    """
    This class is instantiated to represent all products for a route.
    """
    def __init__(self, final_library: Library,
                 template: str = None,
                 hits: str = None,
                 hits_names: List[str] = None,
                 batch_num: int = None,
                 atoms_ids_expansion: dict = None):
        self.products: pd.DataFrame = None
        self.final_library: Library = final_library
        self.output_dir: str = final_library.output_dir
        self.final_products_csv_path: str = None
        # need Fragmenstein information
        self.template: str = template # path to pdb file
        self.hits: str = hits # path to .sdf or .mol file
        self.hits_names: List[str] = hits_names # name of fragments
        self.batch_num: int = batch_num
        self.atoms_ids_expansion: dict = atoms_ids_expansion
        self.placements: pd.DataFrame = None
        self.n_cores: int = 8
        self.timeout: int = 240
        self.output_path: str = None

    def get_products(self):
        """
        Main entry to the Slipper class. This function is used to get the products from the final library.
        """
        slipper_synth = SlipperSynthesizer(self.final_library, self.atoms_ids_expansion)
        self.products: pd.DataFrame = slipper_synth.get_products()
        if self.atoms_ids_expansion is not None:
            slipper_synth.label_products()
        slipper_synth.save_products()
        self.final_products_csv_path: str = slipper_synth.final_products_csv_path
        return self.products

    def place_products(self):
        """
        This function is used to place the products with Fragmenstein.
        """
        slipper_fitter = SlipperFitter(self.products, self.template, self.hits, self.hits_names, self.n_cores,
                                       self.timeout, self.batch_num, self.output_dir, self.final_products_csv_path)
        self.placements: pd.DataFrame = slipper_fitter.fit_products()
        slipper_fitter.save_placements()
        self.output_path: str = slipper_fitter.output_path
        return self.placements

    def _should_delete_file(self, file, suffixes_to_keep):
        """
        Determine if a file should be deleted based on its suffix.
        """
        return not any(file.endswith(suffix) for suffix in suffixes_to_keep)

    def _delete_file_or_directory(self, path):
        """
        Delete a file or directory at the given path.
        """
        try:
            if os.path.isfile(path) or os.path.islink(path):
                os.unlink(path)
                print(f"Deleted file: {path}")
            elif os.path.isdir(path):
                shutil.rmtree(path)
                print(f"Deleted directory: {path}")
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (path, e))

    def clean_up_placements(self):
        """
        This function is used to remove extra files that are generated by Fragmenstein.
        """
        suffixes_to_keep = ['.minimised.json', '.holo_minimised.pdb', '.minimised.mol', '.csv']
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

    def _should_delete_file(self, file, suffixes_to_keep):
        """
        Determine if a file should be deleted based on its suffix.
        """
        return not any(file.endswith(suffix) for suffix in suffixes_to_keep)

    def _delete_file_or_directory(self, path):
        """
        Delete a file or directory at the given path.
        """
        try:
            if os.path.isfile(path) or os.path.islink(path):
                os.unlink(path)
                print(f"Deleted file: {path}")
            elif os.path.isdir(path):
                shutil.rmtree(path)
                print(f"Deleted directory: {path}")
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (path, e))

    def clean_up_placements(self):
        """
        This function is used to remove extra files that are generated by Fragmenstein.
        """
        suffixes_to_keep = ['.minimised.json', '.holo_minimised.pdb', '.minimised.mol', '.csv']
        for root, dirs, files in os.walk(self.output_path):
            for file in files:
                if self._should_delete_file(file, suffixes_to_keep):
                    file_path = os.path.join(root, file)
                    self._delete_file_or_directory(file_path)

