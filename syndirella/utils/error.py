#!/usr/bin/env python3
"""
error.py

Possible errors when using syndirella.
"""
from rdkit import Chem


class ChemicalErrorBase(Exception):
    def __init__(self,
                 message: str,
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None):
        self.message = message
        self.inchi = inchi
        self.smiles = smiles
        self.mol = mol
        super().__init__(self.message)


class MolError(ChemicalErrorBase):
    def __init__(self,
                 route_uuid: str | None = None,
                 message: str = "Could not create a Chem.Mol object.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class SMARTSError(ChemicalErrorBase):
    def __init__(self,
                 route_uuid: str | None = None,
                 message: str = "An error occurred with SMARTS handling.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class NoReactants(ChemicalErrorBase):
    def __init__(self,
                 route_uuid: str | None = None,
                 message: str = "No reactants found for the reaction.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class NoSynthesisRoute(ChemicalErrorBase):
    def __init__(self,
                 message: str = "No synthesis route could be found.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None):
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class ReactionError(ChemicalErrorBase):
    def __init__(self,
                 route_uuid: str,
                 message: str = "An error occurred during the reaction.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class ProductFormationError(ChemicalErrorBase):
    def __init__(self,
                 route_uuid: str,
                 message: str = "Failed to form the scaffold.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class ScaffoldPlacementError(ChemicalErrorBase):
    def __init__(self,
                 route_uuid: str | None = None,
                 message: str = "Scaffold could not be successfully placed.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class PlacementError(ChemicalErrorBase):
    def __init__(self,
                 route_uuid: str,
                 message: str = "Placement could not be successfully performed.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class NoScaffold(ChemicalErrorBase):
    def __init__(self,
                 route_uuid: str,
                 message: str = "No scaffold could be found.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class NoToHippo(ChemicalErrorBase):
    def __init__(self,
                 route_uuid: str,
                 message: str = "The to_hippo pickle file could not be successfully created.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class APIQueryError(ChemicalErrorBase):
    def __init__(self,
                 message: str = "The API query did not successfully return critical information.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None,
                 route_uuid: str | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class APIRetryLimitExceeded(ChemicalErrorBase):
    def __init__(self,
                 message: str = "The API query retry limit has been exceeded.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None,
                 route_uuid: str | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class SingleReactantElabError(ChemicalErrorBase):
    def __init__(self,
                 message: str = "There was an error with producing a series with only a single reactant elaborated.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None,
                 route_uuid: str | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class AiZynthFinderError(ChemicalErrorBase):
    def __init__(self,
                 message: str = "AiZynthFinder could not be successfully run.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None,
                 route_uuid: str | None = None):
        self.route_uuid = route_uuid
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)


class USPTOTemplateValidationError(ChemicalErrorBase):
    def __init__(self,
                 message: str = "USPTO template validation failed.",
                 inchi: str | None = None,
                 smiles: str | None = None,
                 mol: Chem.Mol | None = None,
                 route_uuid: str | None = None,
                 num_routes_found: int = 0,
                 num_routes_with_valid_templates: int = 0):
        self.route_uuid = route_uuid
        self.num_routes_found = num_routes_found
        self.num_routes_with_valid_templates = num_routes_with_valid_templates
        super().__init__(message=message, inchi=inchi, smiles=smiles, mol=mol)
