#!/usr/bin/env python3
"""
Route module containing base configuration classes and route-related objects.
"""

from dataclasses import dataclass
from typing import Optional

from syndirella.constants import DatabaseSearchTool, DEFAULT_DATABASE_SEARCH_TOOL


class BaseRouteConfig:
    """Base configuration class for route objects to reduce variable definitions."""
    def __init__(self, output_dir: str, id: str, atom_diff_min: int, atom_diff_max: int,
                 filter: bool = False, elab_single_reactant: bool = False,
                 elab_single_reactant_int: Optional[int] = None, db_search_tool: DatabaseSearchTool = DEFAULT_DATABASE_SEARCH_TOOL):
        self.output_dir = output_dir
        self.id = id
        self.atom_diff_min = atom_diff_min
        self.atom_diff_max = atom_diff_max
        self.filter = filter
        self.elab_single_reactant = elab_single_reactant
        self.elab_single_reactant_int = elab_single_reactant_int
        self.db_search_tool = db_search_tool
        
        # Validate configuration after initialization
        if self.elab_single_reactant and self.elab_single_reactant_int is None:
            self.elab_single_reactant_int = 0


@dataclass
class LibraryConfig:
    """Configuration class for Library objects."""
    output_dir: str
    id: str
    atom_diff_min: int
    atom_diff_max: int
    num_steps: int
    current_step: int
    route_uuid: str
    filter: bool = False
    elab_single_reactant: bool = False
    elab_single_reactant_int: Optional[int] = None
    db_search_tool: DatabaseSearchTool = DEFAULT_DATABASE_SEARCH_TOOL
    reference_db: str = None


@dataclass
class WorkshopConfig:
    """Configuration class for CobblersWorkshop objects."""
    output_dir: str
    id: str
    atom_diff_min: int
    atom_diff_max: int
    product: str
    reactants: list
    reaction_names: list
    num_steps: int
    filter: bool = False
    elab_single_reactant: bool = False
    elab_single_reactant_int: Optional[int] = None
    db_search_tool: DatabaseSearchTool = DEFAULT_DATABASE_SEARCH_TOOL
    atoms_ids_expansion: Optional[dict] = None


@dataclass
class BenchConfig:
    """Configuration class for CobblerBench objects."""
    output_dir: str
    id: str
    atom_diff_min: int
    atom_diff_max: int
    reaction_name: str
    reactant_smiles: tuple
    filter: bool = False
    elab_single_reactant: bool = False
    elab_single_reactant_int: Optional[int] = None
    db_search_tool: DatabaseSearchTool = DEFAULT_DATABASE_SEARCH_TOOL
