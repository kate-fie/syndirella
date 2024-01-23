#!/usr/bin/env python3
"""
syndirella.cobblers_workshop._base.py

This module contains the CobblersWorkshop class. One instance of this object is used to describe a full route.
"""

from ._cobbler_bench import CobblerBench
from ._database_search import DatabaseSearch
from typing import (List, Dict, Tuple, Union, Optional)
from syndirella.cobblers_workshop._library import Library
from syndirella.smarts import SMARTSHandler
from rdkit import Chem
from rdkit.Chem import rdinchi
from syndirella.slipper._base import Slipper


class CobblersWorkshop():
    """
    This is the CobblersWorkshop class. It represents a full route.
    """
    def __init__(self, product: str, reactants: List[Tuple], reaction_names: List[str],
                 num_steps: int, output_dir: str, database_search: DatabaseSearch):
        self.product: str = product
        self.id: str = self.generate_inchi_ID()
        self.reactants: List[Tuple[str]] = reactants
        self.reaction_names: List[str] = reaction_names
        self.num_steps: int = num_steps
        self.output_dir: str = output_dir
        self.database_search: DatabaseSearch = database_search
        self.smarts_handler = SMARTSHandler()
        self.cobbler_benches: List[CobblerBench] = None # is this actually useful?
        self.first_library: Library = None
        self.final_library: Library = None

    def generate_inchi_ID(self) -> str:
        """
        This function is used to generate a unique id for the route just using the product.
        """
        ID = rdinchi.MolToInchi(Chem.MolFromSmiles(self.product))
        id = rdinchi.InchiToInchiKey(ID[0])
        return id

    def get_final_library(self):
        """
        This function is used to get the final library of products. It is the main function that is called.
        """
        if self.num_steps == 1:
            print(f"There is 1 step in this route using {self.reaction_names[0]}")
            self.get_final_library_one_step()
        if self.num_steps == 2:
            print(f"There are 2 steps in this route using {self.reaction_names[0]} then {self.reaction_names[1]}")
            self.get_final_library_two_steps()
        if self.num_steps > 2:
            raise NotImplementedError("Routes with more than 2 steps are not yet implemented.")
        return self.final_library

    def get_final_library_one_step(self):
        """
        This function is used to get the final library of products for a one step route.
        """
        # Need to convert to a single step
        reactants: Tuple[str] = self.reactants[0]
        reaction_name: str = self.reaction_names[0]
        current_step = 1
        cobbler_bench = CobblerBench(self.product, reactants, reaction_name, self.output_dir, self.smarts_handler,
                                     self.id, self.num_steps, current_step)
        self.final_library = cobbler_bench.find_analogues_first_step()

    def get_final_library_two_steps(self):
        """
        This function is used to get the final library of products for a two step route.
        """
        # Need to convert to a two-step
        reactants1 = self.reactants[0]
        reactants2 = self.reactants[1]
        reaction_name1 = self.reaction_names[0]
        reaction_name2 = self.reaction_names[1]
        current_step = 1
        cobbler_bench1 = CobblerBench(self.product, reactants1, reaction_name1, self.output_dir, self.smarts_handler,
                                      self.id, self.num_steps, current_step)
        self.first_library = cobbler_bench1.find_analogues_first_step()
        first_slipper = Slipper(self.first_library)
        first_slipper.get_products()

        current_step = 2
        cobbler_bench2 = CobblerBench(self.product, reactants2, reaction_name2, self.output_dir, self.smarts_handler,
                                      self.id, self.num_steps, current_step)
        self.final_library = cobbler_bench2.find_analogues_first_step()
