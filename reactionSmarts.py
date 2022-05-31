#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 13:26:19 2021

@author: ruben
"""
import itertools

from collections import OrderedDict

from chemUtils.distances import findCloseAtomIdxsFromRefAtomIdxs
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
    return tuple([(key, Chem.MolFromSmarts(val)) for key,val in smarts_dict.items()])

PATTERN_PRODUCTS= from_SMARTS_to_patterns(_product_smarts)
PATTERN_PRODUCTS_DICT = OrderedDict(PATTERN_PRODUCTS)
PATTERN_REACTANT1= from_SMARTS_to_patterns(_reactant_smarts1)
PATTERN_REACTANT2= from_SMARTS_to_patterns(_reactant_smarts2)


def _checkOneReactionSmartInAttachment(mol, patt, attachment_region_idxs): #This is not enough, it does not guarantee that the molecule was modified at the attachment point
    #TODO: improve this method to check which bonds were broken

    matched_indices = mol.GetSubstructMatches(patt)
    matching = bool(sum([1 if attachment_region_idxs.intersection(option) else 0 for option in matched_indices]))
    return matching

    # matched_indices = mol.GetSubstructMatches(patt)
    # for option in matched_indices:
    #     if len(keepCloseIdxsToAtomIdxs(mol, option, attachment_region_idxs))>0:
    #         return True
    # return False

def checkReactionSmartInAttachment(mol, attachmentIdxs, valid_patterns,
                                       n_hops_from_attachment=config.SMARTS_N_HOPS_FROM_ATTACHMENT):
  if mol is None or attachmentIdxs is None:
      return [ (name, False) for name, patt in valid_patterns]

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



def runReactionInAttachment(reactant1, reactant2, reactionName, refMol, ref_attachmentIdxs,
                            n_hops_from_attachment=config.SMARTS_N_HOPS_FROM_ATTACHMENT):

    ref_attachmentIdxs = expand2DIdxsToNeigs(refMol, ref_attachmentIdxs, n_hops_from_attachment)


    rxn = reactions_rxns[reactionName]
    products = rxn.RunReactants([reactant1, reactant2]) + rxn.RunReactants([reactant2, reactant1]) #itertools.permutations

    if len(products) == 0:
        return []
    products = list(itertools.chain.from_iterable(products))
    list(map(Chem.SanitizeMol, products))
    products = {Chem.MolToSmiles(mol): mol for mol in products if mol is not None} #To remove duplicates
    products = list(products.values())
    products = myMap(lambda molProd: embedMolUsingRefs(molProd, [reactant1, reactant2]), products)

    def checkIfValidReaction(molProd):
        if molProd is None:
            return False
        # Chem.SanitizeMol(molProd)
        # molProd.UpdatePropertyCache()
        #TODO: findCloseAtomsToAtomIdx should be able to find closest to original??
        attachmentIdxs =  [findCloseAtomIdxsFromRefAtomIdxs(molProd, refMol, idx) for idx in ref_attachmentIdxs]
        attachmentIdxs = filter(None.__ne__, attachmentIdxs)
        attachmentIdxs = set(itertools.chain.from_iterable(attachmentIdxs))
        # from molPlot import plotMols; plotMols([molProd, refMol])
        return _checkOneReactionSmartInAttachment(molProd, PATTERN_PRODUCTS_DICT[reactionName], attachmentIdxs)

    isProdValid = myMap(checkIfValidReaction, products)
    return [prod for i, prod in enumerate(products) if isProdValid[i] ]

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