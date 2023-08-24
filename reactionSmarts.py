#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 13:26:19 2021

@author: ruben
"""
import itertools

from collections import OrderedDict

from chemUtils.geometry import findCloseAtomIdxsFromRefAtomIdxs
from rdkit import Chem
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

from config import config
from constants import reaction_smarts
from embed import embedMolUsingRefs
from engine import myMap
from findNeighbours import expand2DIdxsToNeigs

reactions_rxns = {name: ReactionFromSmarts(val) for name, val in reaction_smarts.items()}

_reactant_smarts1 = OrderedDict()
_reactant_smarts2 = OrderedDict()
_product_smarts = OrderedDict()
for name, smart in reaction_smarts.items():
    reactants, prod = smart.split(">>")
    react1, react2 = reactants.split(".")
    _reactant_smarts1[name] = react1
    _reactant_smarts2[name] = react2
    _product_smarts[name] = prod


def from_SMARTS_to_patterns(smarts_dict):
    return tuple([(key, Chem.MolFromSmarts(val)) for key, val in smarts_dict.items()])


PATTERN_PRODUCTS = from_SMARTS_to_patterns(_product_smarts)
PATTERN_PRODUCTS_DICT = OrderedDict(PATTERN_PRODUCTS)
PATTERN_REACTANT1 = from_SMARTS_to_patterns(_reactant_smarts1)
PATTERN_REACTANT2 = from_SMARTS_to_patterns(_reactant_smarts2)

REACTANT1_DICT = _reactant_smarts1
REACTANT2_DICT = _reactant_smarts2

# SPECIFIC REACTION MATCHING
def _contains_ketone_oxygen(mol, atom_indices):
    """
    Check if a list of atom indices contains an index of a ketone oxygen.

    :param mol: An RDKit molecule object.
    :param atom_indices: List of atom indices to check.
    :return: True if ketone oxygen index is found, False otherwise.
    """
    for idx in atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1:  # Check if it's an oxygen with one connection.
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(idx, neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    # Found a carbon double-bonded to the oxygen.
                    return True
    return False

def checkReactionSmartInAllAtomsAndReactants(reactant1, reactant1_attachment_idx, reactant2, reactant2_attachment_idx, reaction_name):
    """
    Given two reactants, their attachment index in product, and reaction name, returns the reaction atoms involved in the reaction.
    Under assumption the reaction name is correct to the reactants.

    :param reactant1:
    :param reactant2:
    :param reaction_name:
    :return:
    """
    # Check if pattern of reactant 1 is in one of the reactant molecules at the attachment point, if not, check reactant 2
    reactant1_attachment_idx = set(reactant1_attachment_idx)
    reactant2_attachment_idx = set(reactant2_attachment_idx)
    matched_atoms = {"reactant1": {}, "reactant2": {}} # tuple of (reactant_smarts, mol, matched_indices)

    reactant_smarts = {"reactant1": REACTANT1_DICT[reaction_name], "reactant2": REACTANT2_DICT[reaction_name]}

    # check reactant1 input exhausitvely against both reactant smart patterns
    patt1 = Chem.MolFromSmarts(reactant_smarts["reactant1"])
    patt2 = Chem.MolFromSmarts(reactant_smarts["reactant2"])
    matches1 = reactant1.GetSubstructMatches(patt1)
    matches2 = reactant1.GetSubstructMatches(patt2)
    found = False
    if matches1:
        for i in range(len(matches1)):
            matching = bool(sum([1 if reactant1_attachment_idx.intersection(option) else 0 for option in matches1]))
            if matching: # SHOULD BE EXACTLY MATCHING SINCE LOOKING AT ATTACHMENT INDEX
                matched_atoms["reactant1"] = (reactant_smarts["reactant1"], reactant1, matches1[i])
                found = True
    if found is False and matches2:
        for i in range(len(matches2)):
            matching = bool(sum([1 if reactant1_attachment_idx.intersection(option) else 0 for option in matches2]))
            if matching: # SHOULD BE EXACTLY MATCHING SINCE LOOKING AT ATTACHMENT INDEX
                matched_atoms["reactant1"] = (reactant_smarts["reactant2"], reactant1, matches2[i])

    # check reactant2 input exhausitvely against both reactant smart patterns
    matches1 = reactant2.GetSubstructMatches(patt1)
    matches2 = reactant2.GetSubstructMatches(patt2)
    found = False
    if matches1:
        for i in range(len(matches1)):
            matching = bool(sum([1 if reactant2_attachment_idx.intersection(option) else 0 for option in matches1]))
            if matching:
                matched_atoms["reactant2"] = (reactant_smarts["reactant1"], reactant2, matches1[i])
                found = True
    if found is False and matches2:
        for i in range(len(matches2)):
            matching = bool(sum([1 if reactant2_attachment_idx.intersection(option) else 0 for option in matches2]))
            if matching:
                matched_atoms["reactant2"] = (reactant_smarts["reactant2"], reactant2, matches2[i])

    # check that there is at least one item in the matched_atoms dict
    if len(matched_atoms["reactant1"]) == 0 or len(matched_atoms["reactant2"]) == 0:
        print('WARNING: No atoms found involved in reaction ' + reaction_name + ' in mol ' + Chem.MolToSmiles(reactant1) + ' and ' + Chem.MolToSmiles(reactant2))
        return None

    return matched_atoms

def checkSpecificReactionSmartInReactant(smiles, reaction_name, reaction_smarts):
    """
    Given a molecule and a reaction smarts that it should contain because it is a superstructure of it, checks if the molecule has the reaction pattern.

    :param mol:
    :param reaction_pattern:
    :return:
    """
    mol = Chem.MolFromSmiles(smiles)
    Chem.AddHs(mol)
    matched_reactions = []
    matching = _checkOneReactionSmartInAttachmentSTRICT_noidx(mol, reaction_name, reaction_smarts)
    matched_reactions.append((reaction_name, matching))
    return matched_reactions

def _checkOneReactionSmartInAttachmentSTRICT_noidx(mol, name, smarts):
    patt = Chem.MolFromSmarts(smarts)

    if name == "Amidation":
        matched_indices = mol.GetSubstructMatches(patt)
        if matched_indices:
            if patt.GetNumAtoms() == 3: # carboxylic acid
                # check length == 3
                for i in range(len(matched_indices)): # could have multiple matches ...
                    # check you have hydroxy group in acid
                    if len(matched_indices[i]) == 3:
                        return True
            else: # primary or secondary amine
                # check it is not an amide, so does not contain an oxygen
                for i in range(len(matched_indices)): # could have multiple matches ...
                    if len(matched_indices[i])==1:
                        return True

    if name == "Amide_schotten-baumann":
        pass
    if name == "Reductive_amination":
        pass
    if name == "N-nucleophilic_aromatic_substitution":
        pass
    if name == "Sp2-sp2_Suzuki_coupling":
        pass
    if name == "Formation_of_urea_from_two_amines":
        pass

    return False

def _checkOneReactionSmartInAttachment(mol, patt,
                                       attachment_region_idxs):  # This is not enough, it does not guarantee that the molecule was modified at the attachment point
    # TODO: improve this method to check which bonds were broken

    matched_indices = mol.GetSubstructMatches(patt)
    matching = bool(sum([1 if attachment_region_idxs.intersection(option) else 0 for option in matched_indices]))
    return matching

    # matched_indices = mol.GetSubstructMatches(patt)
    # for option in matched_indices:
    #     if len(keepCloseIdxsToAtomIdxs(mol, option, attachment_region_idxs))>0:
    #         return True
    # return False

def _checkOneReactionSmartInAttachmentSTRICT(mol, name, patt, attachment_region_idxs, reactantNum):
    """
    Strict reaction smart check, specifics of the reaction are checked. Cons: Manual.

    :param mol:
    :param patt: either reactant1 or reactant2 pattern
    :param attachment_region_idxs:
    :return:
    """
    if name == "Amidation":
        matched_indices = mol.GetSubstructMatches(patt)
        matching = bool(sum([1 if attachment_region_idxs.intersection(option) else 0 for option in matched_indices]))
        if matching:
            if reactantNum == 1: # carboxylic acid
                # check length == 3
                if len(matched_indices) != 3:
                    return False
            else: # primary or secondary amine
                # check it is not an amide, so does not contain an oxygen
                amide = _contains_ketone_oxygen(mol, matching)
                if amide:
                    return False

    if name == "Amide_schotten-baumann":
        pass
    if name == "Reductive_amination":
        pass
    if name == "N-nucleophilic_aromatic_substitution":
        pass
    if name == "Sp2-sp2_Suzuki_coupling":
        pass
    if name == "Formation_of_urea_from_two_amines":
        pass

    return True


def checkReactionSmartAroundAttachment(mol, attachmentIdxs, valid_patterns, reactantNum,
                                       n_hops_from_attachment=config.SMARTS_N_HOPS_FROM_ATTACHMENT):
    """
    Given a molecule and an attachment point, checks if the molecule has a reaction pattern around the attachment point.
    Implementing strict check, this is where streochemistry and electronics can be checked.

    :param mol_attachmentIdx:
    :param valid_patterns:
    :param n_hops_from_attachment:
    :return:
    """
    if mol is None or attachmentIdxs is None:
        return [(name, False) for name, patt in valid_patterns]

    # I don't think we need to expand
    attachment_region_idxs = expand2DIdxsToNeigs(mol, attachmentIdxs, n_hops_from_attachment)

    matched_reactions = []
    for name, patt in valid_patterns:
        matching = _checkOneReactionSmartInAttachmentSTRICT(mol, name, patt, attachment_region_idxs, reactantNum)
        matched_reactions.append((name, matching))

    return matched_reactions

def checkReactionSmartInAttachment(mol, attachmentIdxs, valid_patterns,
                                   n_hops_from_attachment=config.SMARTS_N_HOPS_FROM_ATTACHMENT):
    if mol is None or attachmentIdxs is None:
        return [(name, False) for name, patt in valid_patterns]

    attachment_region_idxs = expand2DIdxsToNeigs(mol, attachmentIdxs, n_hops_from_attachment)

    # print(attachment_region_idxs)
    # from molPlot import plotMols
    # plotMols([mol])
    matched_reactions = []
    for name, patt in valid_patterns:
        matching = _checkOneReactionSmartInAttachment(mol, patt, attachment_region_idxs)
        matched_reactions.append((name, matching))
        # print(name, matched_indices, matching)
    # from molPlot import plotMols
    # plotMols([mol])
    return matched_reactions


def checkReactionSmartInAllAtoms(mol, reaction_name) -> list:

    matched_reactions = {"reactant1": {}, "reactant2": {}}

    PATTERNS = {"reactant1": PATTERN_REACTANT1, "reactant2": PATTERN_REACTANT2}

    for reactant, patterns in PATTERNS.items():
        for name, patt in patterns:
            matches = mol.GetSubstructMatches(patt)
            if matches:  # If there are matches
                matched_reactions[reactant][name] = (True, matches)
            else:  # If there are no matches
                matched_reactions[reactant][name] = (False, [])
    if reaction_name is not None:
        try:
            matched_reactions1 = matched_reactions['reactant1'][reaction_name]
            matched_reactions2 = matched_reactions['reactant2'][reaction_name]
            if matched_reactions1[0]:  # If reaction1 has matches
                return matched_reactions1[1]
            elif matched_reactions2[0]:  # If reaction2 has matches
                return matched_reactions2[1]
            else:  # If neither reaction has matches
                print('WARNING: No atoms found involved in reaction ' + reaction_name + ' in mol ' + Chem.MolToSmiles(mol))
                return []
        except KeyError:
            print('WARNING: Reaction name not in smarts dictionary.')
            return None
    else:
        print('WARNING: Reaction name not specified.')
        return None


def runReactionInAttachment(reactant1, reactant2, reactionName, refMol, ref_attachmentIdxs,
                            n_hops_from_attachment=config.SMARTS_N_HOPS_FROM_ATTACHMENT):
    ref_attachmentIdxs = expand2DIdxsToNeigs(refMol, ref_attachmentIdxs, n_hops_from_attachment)

    rxn = reactions_rxns[reactionName]
    products = rxn.RunReactants([reactant1, reactant2]) + rxn.RunReactants(
        [reactant2, reactant1])  # itertools.permutations

    if len(products) == 0:
        return []
    products = list(itertools.chain.from_iterable(products))
    list(map(Chem.SanitizeMol, products))
    products = {Chem.MolToSmiles(mol): mol for mol in products if mol is not None}  # To remove duplicates
    products = list(products.values())
    products = myMap(lambda molProd: embedMolUsingRefs(molProd, [reactant1, reactant2]), products)

    def checkIfValidReaction(molProd):
        if molProd is None:
            return False
        # Chem.SanitizeMol(molProd)
        # molProd.UpdatePropertyCache()
        # TODO: findCloseAtomsToAtomIdx should be able to find closest to original??
        attachmentIdxs = [findCloseAtomIdxsFromRefAtomIdxs(molProd, refMol, idx) for idx in ref_attachmentIdxs]
        attachmentIdxs = filter(None.__ne__, attachmentIdxs)
        attachmentIdxs = set(itertools.chain.from_iterable(attachmentIdxs))
        # from molPlot import plotMols; plotMols([molProd, refMol])
        return _checkOneReactionSmartInAttachment(molProd, PATTERN_PRODUCTS_DICT[reactionName], attachmentIdxs)

    isProdValid = myMap(checkIfValidReaction, products)
    return [prod for i, prod in enumerate(products) if isProdValid[i]]

# def runReactionInAttachment_v1(reactant1, reactant2, reactionName, attachmentIdxs1,  attachmentIdxs2, referenceMol,
#                             n_hops_from_attachment=SMARTS_N_HOPS_FROM_ATTACHMENT):
#
#     attachment_region_idxs1 = expand2DIdxsToNeigs(reactant1, attachmentIdxs1, n_hops_from_attachment)
#     attachment_region_idxs2 = expand2DIdxsToNeigs(reactant2, attachmentIdxs2, n_hops_from_attachment)
#     rxn = reactions_rxns[reactionName]
#     products = rxn.RunReactants([reactant1, reactant2])
#
#     # if len(products) == 0:
#     #     return None
#     products = myMap(lambda molProd: embedMolByMCS(molProd, referenceMol), products[0])
#
#     def checkIfValidReaction(molProd):
#         if molProd is None:
#             return False
#         Chem.SanitizeMol(molProd)
#         molProd.UpdatePropertyCache()
#         #TODO: findCloseAtomsToAtomIdx should be able to find closest to original??
#         attachmentIdxs =  [findCloseAtomIdxsFromRefAtomIdx(molProd, reactant1, idx) for idx in attachment_region_idxs1]
#         attachmentIdxs += [findCloseAtomIdxsFromRefAtomIdx(molProd, reactant2, idx) for idx in attachment_region_idxs2]
#         # attachmentIdxs =  [findCloseAtomIdxsToAtomIdx(referenceMol, molProd, idx) for idx in attachment_region_idxs1]
#         # attachmentIdxs += [findCloseAtomIdxsToAtomIdx(referenceMol, molProd, idx) for idx in attachment_region_idxs2]
#         attachmentIdxs = filter(None.__ne__, attachmentIdxs)
#         attachmentIdxs = set(itertools.chain.from_iterable(attachmentIdxs))
#         return _checkOneReactionSmartInAttachment(molProd, PATTERN_PRODUCTS_DICT[reactionName], attachmentIdxs)
#
#     isProdValid = myMap(checkIfValidReaction, products)
#     return [prod for i, prod in enumerate(products) if isProdValid[i] ]
