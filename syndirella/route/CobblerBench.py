#!venv/bin/env python3
"""
syndirella.CobblerBench.py

This module contains the CobblerBench class. One instance is made for each step.
"""

import logging
from typing import (List, Tuple)

from rdkit import Chem

from .Library import Library
from .Reaction import Reaction
from .SMARTSHandler import SMARTSHandler
from syndirella.utils.error import ReactionError, SMARTSError
from ..constants import DatabaseSearchTool, DEFAULT_DATABASE_SEARCH_TOOL


class CobblerBench:
    """
    This class is used to perform the whole process of finding analogues of reactants. Since the final elaborated products
    are 'slippers' in this analogy, the CobblerBench is where these slippers are made.

    It is given a step of a full route.
    """

    def __init__(self,
                 scaffold: str,
                 reaction_name: str,
                 reactant_smiles: Tuple[str] | str,
                 output_dir: str,
                 filter: bool,
                 id: str,
                 num_steps: int,
                 current_step: int,
                 atom_diff_min: int,
                 atom_diff_max: int,
                 elab_single_reactant: bool,
                 route_uuid: str,
                 elab_single_reactant_int: int | None = None,
                 db_search_tool: DatabaseSearchTool = DEFAULT_DATABASE_SEARCH_TOOL,
                 reference_db: str | None = None):
        self.product: Chem.Mol = Chem.MolFromSmiles(scaffold)
        self.reaction_name: str = reaction_name
        self.reactant_smiles: Tuple[str] | str = reactant_smiles
        self.reactants: Tuple[Chem.Mol] = self._make_reactant_mols(reactant_smiles)
        self.output_dir: str = output_dir
        self.smarts_handler: SMARTSHandler = SMARTSHandler()
        self.id: str = id
        self.filter: bool = filter
        self.atom_diff_min: int = atom_diff_min
        self.atom_diff_max: int = atom_diff_max
        self.elab_single_reactant: bool = elab_single_reactant
        self.elab_single_reactant_int: int | None = elab_single_reactant_int
        self.db_search_tool: DatabaseSearchTool = db_search_tool
        self.num_steps: int = num_steps
        self.current_step: int = current_step
        self.route_uuid: str = route_uuid
        self.reference_db: str | None = reference_db
        
        self.reaction: Reaction = None
        self.library = None
        self.additional_reactions: List[Tuple[str, Tuple[str, str]]] = ()
        self.logger = logging.getLogger(f"{__name__}")

        self.define_reaction()  # This is where the Reaction object is created

    def _make_reactant_mols(self, reactants) -> Tuple[Chem.Mol]:
        """
        This function is used to make reactant molecules from the input SMILES.
        """
        if type(reactants) == str:
            reactant_mols = [Chem.MolFromSmiles(reactants)]
        else:
            reactant_mols = [Chem.MolFromSmiles(reactant) for reactant in reactants]
        return tuple(reactant_mols)

    def check_reaction(self):
        """
        This function is used to check the reaction is valid.
        """
        if self.reaction_name not in self.smarts_handler.reaction_smarts:
            self.logger.error(f"Reaction name {self.reaction_name} not found in SMARTS handler.")
            raise SMARTSError(message=f"Reaction name {self.reaction_name} not found in SMARTS handler.",
                              smiles=self.product,
                              route_uuid=self.route_uuid,
                              inchi=self.id)
        reaction = Reaction(scaffold=self.product,
                            reactants=self.reactants,
                            reaction_name=self.reaction_name,
                            smarts_handler=self.smarts_handler,
                            route_uuid=self.route_uuid)
        return reaction

    def define_reaction(self):
        """
        This function is used to define the reaction and finding ids and attachment atoms.
        """
        try:
            self.reaction = self.check_reaction()
        except (ValueError, SMARTSError) as e:
            self.logger.error(f"Reaction {self.reaction_name} could not be defined.")
            if isinstance(e, ValueError):
                message = e.args[0]
            elif isinstance(e, SMARTSError):
                message = e.message
            else:
                message = f"Reaction {self.reaction_name} could not be defined."
            raise ReactionError(message=message,
                                smiles=self.product,
                                route_uuid=self.route_uuid,
                                inchi=self.id)

    def get_additional_reactions(self) -> bool:
        """
        Gets the additional reaction by outputting edited reactants via SMARTS and new reaction names all controlled
        by the Reaction class.
        """
        new_reactions: List[Tuple[str, Tuple[str, str]]] = self.reaction.get_additional_reactions()
        if len(new_reactions) == 0:
            self.logger.info(f"No additional reactions found for {self.reaction_name}. Continuing with original.")
            return False
        self.additional_reactions = new_reactions
        return True

    def find_analogues(self):
        """
        This function is used to find analogues of any step along the route.
        """
        # Find the analogues of reactants
        self.library = Library(
            reaction=self.reaction,
            output_dir=self.output_dir,
            id=self.id,
            num_steps=self.num_steps,  # Single step for this bench
            current_step=self.current_step,
            filter=self.filter,
            route_uuid=self.route_uuid,  # Will be set by parent
            atom_diff_min=self.atom_diff_min,
            atom_diff_max=self.atom_diff_max,
            db_search_tool=self.db_search_tool,
            elab_single_reactant_int=self.elab_single_reactant_int,
            elab_single_reactant=self.elab_single_reactant,
            reference_db=self.reference_db,
        )
        self.library.create_library()
        return self.library
