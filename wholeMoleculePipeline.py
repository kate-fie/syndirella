#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:12:49 2021

@author: ruben
"""
import itertools
import os
from collections import OrderedDict
from typing import Tuple, Iterable

import numpy as np
import pandas as pd
from chemUtils.geometry.molDistances import findCloseAtomIdxsFromRefAtomIdxs, find_exact_match_atom_indices, find_substruct_atom_match
from chemUtils.substructure import splitExtendedMolBySubstruct
from chemUtils.substructure import find_attachmentIdxs_fromMolAndSubstruct
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.PandasTools import WriteSDF
from rdkit.rdBase import BlockLogs

from config import config
from constants import fromReactionFullNameToReactantName, REACTIONS_NAMES, REACTION_ATOM_NUMS
from databaseSearch.databaseSearch import perform_database_search
from embed import embedMolUsingRefs
from engine import myMap
from reactionSmarts import (checkReactionSmartInAttachment, PATTERN_REACTANT1, PATTERN_REACTANT2, \
    runReactionInAttachment, checkReactionSmartInAllAtoms, reaction_smarts,
                            checkReactionSmartInAllAtomsAndReactants, checkReactionSmartAroundAttachment, checkSpecificReactionSmartInReactant)
from scorer import scorePairOfMols

from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdDistGeom
from rdkit.Chem.rdfiltercatalog import FilterCatalog, FilterCatalogParams

import pickle
from itertools import chain

from PIL import Image, ImageDraw, ImageFont

def pains_filter(df):
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    def is_pain(smile):
        if not isinstance(smile, str):
            return None
        mol = Chem.MolFromSmiles(smile)
        if mol and catalog.HasMatch(mol):
            return True
        return False

    df['is_pains'] = df['product_smi'].apply(is_pain)

    return df

def save_molecule_with_indices(mol, filename, text):
    # Generate a new molecule that includes atom indices
    mol_with_indices = Chem.Mol(mol)
    for atom in mol_with_indices.GetAtoms():
        atom.SetProp('molAtomMapNumber', str(atom.GetIdx()))

    # Convert the molecule with indices to an image
    img = Draw.MolToImage(mol_with_indices, size=(500, 500))

    # Convert RDKit Image to PIL Image
    pil_img = Image.new("RGB", img.size)
    pil_img.paste(img)

    # If there's text to add
    if text:
        draw = ImageDraw.Draw(pil_img)
        draw.text((10, 10), text, (255, 255, 255))  # Coordinates, text, color, and font

    # Save the image to the file
    pil_img.save(filename)

def getFunction_checkStructuralMatch(originalMol):
    originalMol = Chem.AddHs(originalMol, addCoords=True)

    def analyzeOneSmi(simMol):
        # if not simMol:
        #     return None
        logsManager = BlockLogs()
        # Changing mol to explicit Hs which should help with embedding
        simMol = Chem.AddHs(simMol, addCoords=True)
        simEmbedded = embedMolUsingRefs(simMol, originalMol)  # TODO: fix to not return None all the time
        AllChem.EmbedMolecule(originalMol) # Let_s try to embed the original mol, it worked
        # calculate the mw of original and new mol
        mw_original = Descriptors.MolWt(originalMol)
        mw_sim = Descriptors.MolWt(simMol)
        mw_diff = abs(mw_original - mw_sim)
        del logsManager
        if simEmbedded is None:
            return None, None, None, None
        else:
            return scorePairOfMols(originalMol, simEmbedded)[0], simEmbedded, mw_diff, mw_sim

    return analyzeOneSmi


def embedAndScoreSimilarSmiles(refMol: Chem.Mol, smilesList) -> Tuple[Iterable[float], Iterable[Chem.Mol]]:
    """

    :param refMol: The reference mol to embed the smiles
    :param smilesList: A list of smiles
    :return: a list of tuples (structScore, smiles, similarity)
    """

    similar_mols = myMap(Chem.MolFromSmiles, smilesList)
    results = myMap(getFunction_checkStructuralMatch(refMol), similar_mols)

    structuralScore, embedded_similar_mols, mw_diff, mw = zip(*results) if results else ([], [], [], [])
    return structuralScore, embedded_similar_mols, mw_diff, mw


def matchMolToReactan(df, reaction_name, reaction_smarts, reactant_Num):
    """
    Given a dataframe with a column "mol" and a column "attachmentIdx", this function checks if the molecule is a
    reactant of specific reaction provided with reaction_name.

    :param df:
    :param reactantNum:
    :return:
    """
    # if reactantNum == 1:
    #     valid_patterns = PATTERN_REACTANT1
    # elif reactantNum == 2:
    #     valid_patterns = PATTERN_REACTANT2  # Need to have two reactants, flip the order to find both combinations
    # else:
    #     raise ValueError("Only 1 and 2 valid for reactant_num")
    # matching_reactions = myMap(lambda mol_attachmentIdx:
    #                            checkReactionSmartAroundAttachment(*mol_attachmentIdx,
    #                                                           valid_patterns=valid_patterns, reactantNum=reactantNum),
    #                            zip(df["mol"], df["attachmentIdx"]))
    matching_reactions = myMap(lambda smiles:
                               checkSpecificReactionSmartInReactant(smiles, reaction_name, reaction_smarts),
                               df["smi"])

    selected_reactions = {name: [] for name in REACTIONS_NAMES}
    # can return false here
    for i, reaction_name_flag in enumerate(matching_reactions):
        for name, flag in reaction_name_flag:
            selected_reactions[name].append(flag)

    vals = selected_reactions[reaction_name]
    assert len(vals) == len(df), ("Error, the number of values in the selected reactions is not the same as the number "
                                  "of molecules")
    df[fromReactionFullNameToReactantName(name, reactant_Num)] = vals  # making new column in dataframe

    return df


def findReactionAtoms(ori_mol, reactant, reaction_name, results_dir, react_num, attachmentIdxs):
    """
    Given a molecule and a reactant, this function finds the atoms that are involved in the reaction for the reactant.
    :param ori_mol:
    :param reactant:
    :param reaction_name:
    :return:
    """
    # TODO: Need to save if first reactant has already been found or not.
    # TODO: NEED TO CLEAN THIS UP AND MAKE IT LOOK AND OPERATE PRETTIER
    matched_reaction_atoms = checkReactionSmartInAllAtoms(reactant, reaction_name, )
    if matched_reaction_atoms is None:
        print('No reaction atoms found for reactant', Chem.MolToSmiles(reactant), 'in reaction', reaction_name)
        return None
    matched_reaction_atoms = tuple(x for sub in matched_reaction_atoms for x in sub)  # This flattens the tuple
    save_molecule_with_indices(reactant, os.path.join(results_dir,
                                                      f"{reaction_name}_reactant{react_num}.png"), matched_reaction_atoms)

    # Check reaction atoms are in attachment idx
    if not set(attachmentIdxs).issubset(set(matched_reaction_atoms)): # attachmentIdxs is always smaller than matched_reaction_atoms (should be?)
        print('Reaction atoms are not in attachment idxs')
        return None
    else:
        # Intersection of the attachmentIdxs and the matched_reaction_atoms, don't neccesarily want this
        if reaction_name == 'Amidation' and len(matched_reaction_atoms) == 3:
            return matched_reaction_atoms
        matched_reaction_atoms = list(set(attachmentIdxs).intersection(set(matched_reaction_atoms)))

    # Check number of matched_reaction_atoms to correspond with correct number in reaction
    num_rxn_atoms = REACTION_ATOM_NUMS[reaction_name]

    match = False
    for i in range(len(num_rxn_atoms)):
        if len(matched_reaction_atoms) == num_rxn_atoms[i]:
            match = True

    if match is False:
        print('Number of reaction atoms does not match the number in the reaction')
        return None
    return matched_reaction_atoms

def findReactionAtoms_bothReactants(originalMol, reactant1, reactant2, ori_reaction, resultsDir, reactant1_attachmentIdxs, reactant2_attachmentIdxs):
    """
    Given a product, two reactants and a labeled reaction between them, need to find the atoms involved in the reaction on the reactants.

    :param originalMol:
    :param reactant1:
    :param reactant2:
    :param ori_reaction:
    :param resultsDir:
    :param reactant1_attachmentIdxs:
    :param reactant2_attachmentIdxs:
    :return:
    """

    # matched_reaction_atoms is a dictionary of 'reactant1' and 'reactant2' with values as atom indices involved in reaction
    matched_reaction_dict = checkReactionSmartInAllAtomsAndReactants(reactant1, reactant1_attachmentIdxs, reactant2, reactant2_attachmentIdxs, reaction_name=ori_reaction)
    # if len(matched_reaction_atoms[0]) == 0 or len(matched_reaction_atoms[1]) == 0:
    #     print('Reaction atoms not found for reactants', Chem.MolToSmiles(reactant1), Chem.MolToSmiles(reactant2), 'in reaction', ori_reaction)
    #     return None
    matched_reaction_atoms1 = matched_reaction_dict['reactant1'][2]
    save_molecule_with_indices(reactant1, os.path.join(resultsDir,
                                                      f"{ori_reaction}_reactant{1}.png"), matched_reaction_atoms1)
    matched_reaction_atoms2 = matched_reaction_dict['reactant2'][2]
    save_molecule_with_indices(reactant2, os.path.join(resultsDir,
                                                      f"{ori_reaction}_reactant{2}.png"), matched_reaction_atoms2)
    matched_reaction_info = ((matched_reaction_dict['reactant1'][0], matched_reaction_dict['reactant1'][2]), (matched_reaction_dict['reactant2'][0], matched_reaction_dict['reactant2'][2]))
    return matched_reaction_info


def createAnaloguesDf(smisDict, structuralScores, mols, mw_diff, mw, similarity, reactant):
    if similarity:
        similarities, metadatas = zip(*smisDict.values())
        result = OrderedDict(structuralScore=structuralScores, mw_diff=mw_diff, mw=mw, similarities=similarities, smi=list(smisDict.keys()),
                             metadata=metadatas, mol=mols)
    else:
        metadatas = list(smisDict.values())
        result = OrderedDict(structuralScore=structuralScores, mw_diff=mw_diff, mw=mw, smi=list(smisDict.keys()),
                             metadata=metadatas, mol=mols)
    result = pd.DataFrame(result)
    result.dropna(axis=0, inplace=True)

    # add original mol to dataframe
    result.loc[-1] = [0, 0, Chem.Descriptors.ExactMolWt(reactant), Chem.MolToSmiles(reactant), None, None]  # adding a row
    result.index = result.index + 1  # shifting index
    result.sort_index(inplace=True)
    result.sort_values(by="structuralScore", inplace=True)
    result.reset_index(drop=True, inplace=True)
    return result


def processSubMolReactant(reactant, reaction_name, reaction_smarts, reaction_atoms, similarityThr, structuralThr, reactant_Num, resultsDir, substructure=True):
    """
    Given a reactant and the atoms involved in the reaction, this function searches for analogues of the reactant
    with Postera's superstructure search.

    :param reactant: mol object of reactant
    :param reaction_atoms: atoms involved in the reaction
    :param similarityThr:
    :param structuralThr:
    :param resultsDir:
    :param substructure: True if searching for superstructure
    :return: DataFrame with analogues of the reactant.
    """
    print('This is the reactant:', Chem.MolToSmiles(reactant))
    print('These are the atoms involved in the reaction:', reaction_atoms)

    # we know which reactant it is
    smisDict = perform_database_search(reactant, substructure=substructure, thr=similarityThr) # substructure is the search for superstructure

    # Save smisDict to a file, to show structures easily
    out_dict = os.path.join(resultsDir, 'smisDict.pickle')
    with open(out_dict, 'wb') as f:
        pickle.dump(smisDict, f)

    structuralScores, similarMols, mw_diff, mw = embedAndScoreSimilarSmiles(reactant, smisDict.keys())

    similarity = True
    if substructure is True: similarity = False

    df = createAnaloguesDf(smisDict, structuralScores, similarMols, mw_diff, mw, similarity, reactant)
    df.to_csv(os.path.join(resultsDir, f'analogues_preprocess_{Chem.MolToSmiles(reactant)}.csv'))

    # Search through elaborated compounds for the reaction atoms as substructure search
    #equiv_attachmentIdxs = myMap(lambda mol: find_substruct_atom_match(mol, reactant, reaction_atoms), df["mol"])

    # Have to create conformers for the elaborated mols and reactants
    # Do you need a 3D conformer? I'll make a 3D one. I think the 3D conformer was already calcuated when they were embedded to compare similarity


    # equiv_attachmentIdxs2 =  myMap(lambda mol: find_exact_match_atom_indices(mol, reactant, reaction_atoms, test=False), df['mol'])
    #
    # # The equiv_attachementIdx2 is
    #
    # df.loc[:, "attachmentIdx"] = equiv_attachmentIdxs2
    # df = df[df["attachmentIdx"].notna()]
    print('This is the length of the dataframe after substructure search:', len(df))
    if structuralThr is not None:
        df = df[df["structuralScore"] > structuralThr]
    df = matchMolToReactan(df, reaction_name, reaction_smarts, reactant_Num) # should be specific reactant
    df = df.reset_index(drop=True)
    df.to_csv(os.path.join(resultsDir, f'analogues_substruct_{Chem.MolToSmiles(reactant)}_{len(df)}.csv'))

    return df

import pandas as pd
from itertools import product

def apply_reaction(row, reaction_name) -> list():
    """
    Applies a reaction to a row of a dataframe. Need to make sure that the reaction SMARTS is flipped so both sides of
    the reaction are tried out. But need to make sure that the reactants are in the right order like how they were originally discovered.
    Since functional groups decorating a ring do not want to be reacted together even if they match the smarts pattern.

    :param row:
    :param reaction_name:
    :return: list of product SMILES
    """
    #print('This is reactant1:', row['reactant1_smi'])
    #print('This is reactant2:', row['reactant2_smi'])

    # Reactant1 is just the first reactant in the SMARTS, and reactant2 is the second reactant in the SMARTS
    smarts_pattern = reaction_smarts[reaction_name]
    reaction = AllChem.ReactionFromSmarts(smarts_pattern)

    reactant1 = Chem.MolFromSmiles(row['reactant1_smi'])
    reactant2 = Chem.MolFromSmiles(row['reactant2_smi'])
    #reactant1 = Chem.MolFromSmiles('CCOC(=O)/C=C1\CCC(Cc2ccccc2)(C(=O)O)C1')
    #reactant2 = Chem.MolFromSmiles('CN(C)C1CCN(C/C(N)=N/O)CC1')
    # TODO: WHY ARE
    products = reaction.RunReactants((reactant1, reactant2))

    # Convert the products to SMILES strings
    product_smiles = [Chem.MolToSmiles(product[0]) for product in products]

    #refMol = Chem.MolFromSmiles('CC1CCC(=CC(=O)NCCN2CCCCC2)C1')
    #ref_attachmentIdxs = [8]
    # TODO: Test if Ruben's way of finding products works better.
    # products = runReactionInAttachment(reactant1, reactant2, reaction_name, refMol, ref_attachmentIdxs)

    # if len(product_smiles) > 1:
    #     print('This is reactant1:', row['reactant1_smi'])
    #     print('This is reactant2:', row['reactant2_smi'])
    #     print('There is more than one product:', product_smiles)
    return product_smiles

def findReactionCombinations(ori_reaction, analogs_reactant1, analogs_reactant2, resultsDir):
    """
    Given a reaction and the analogues of the reactants, this function finds the reaction combinations that are possible.

    :param ori_reaction:
    :param analogs_reactant1:
    :param analogs_reactant2:
    :param resultsDir:
    :return: DataFrame with reactant combinations
    """
    i = ori_reaction

    # Filter both DataFrames based on the react1_ReactionName and react2_ReactionName columns
    filtered_df1 = analogs_reactant1[analogs_reactant1[fromReactionFullNameToReactantName(i, 1)] == True]
    filtered_df2 = analogs_reactant2[analogs_reactant2[fromReactionFullNameToReactantName(i, 2)] == True]

    # Switching reactant numbers
    if filtered_df2.empty and filtered_df1.empty:
        filtered_df1 = analogs_reactant1[analogs_reactant1[fromReactionFullNameToReactantName(i, 2)] == True]
        filtered_df2 = analogs_reactant2[analogs_reactant2[fromReactionFullNameToReactantName(i, 1)] == True]

    # Get the Cartesian product of the filtered DataFrame based on index
    reactant_combinations = pd.MultiIndex.from_product([filtered_df1.index, filtered_df2.index], names=['reactant1', 'reactant2']).to_frame(index=False)

    filtered_df1 = filtered_df1[['structuralScore', 'mw', 'smi', 'metadata']]
    filtered_df2 = filtered_df2[['structuralScore', 'mw', 'smi', 'metadata']]

    # Merge the information based on the indices
    result_df = (
        reactant_combinations
        .merge(filtered_df1, left_on='reactant1', right_index=True)
        .merge(filtered_df2, left_on='reactant2', right_index=True, suffixes=('_reactant1', '_reactant2'))
    )
    result_df = result_df.drop(columns=['reactant1', 'reactant2'])

    # Rename the columns to the desired format
    result_df.rename(columns={
        'smi_reactant1': 'reactant1_smi',
        'structuralScore_reactant1': 'reactant1_structuralScore',
        'mw_reactant1': 'reactant1_mw',
        'metadata_reactant1': 'reactant1_metadata',
        'smi_reactant2': 'reactant2_smi',
        'structuralScore_reactant2': 'reactant2_structuralScore',
        'mw_reactant2': 'reactant2_mw',
        'metadata_reactant2': 'reactant2_metadata',
    }, inplace=True)

    # Add a column of the reacted together SMILES
    # WHY DO I KEEP GETTING THIS ERROR OF MULTIPLE COLUMNS TO SINGLE COLUMN????
    # SOLN: Happens when you literally have no products.
    result_df['product_smi'] = result_df.apply(apply_reaction, reaction_name=i, axis=1)
    # Use the explode method to transform lists in the 'product_smi' column into separate rows
    result_df = result_df.explode('product_smi').reset_index(drop=True)

    # Extracting the product_smi column and flattening the lists
    all_products = result_df['product_smi']

    # Remove repeated rows
    result_df2 = result_df.copy()
    result_df2.drop_duplicates(subset=['reactant1_smi', 'reactant2_smi', 'product_smi'])
    print(f'{len(result_df)-len(result_df2)} product duplicates removed.')

    # How many unique products are there?
    print(f'{result_df2["product_smi"].nunique()} unique products found.')

    # ADDING PAINS FILTER
    pains_filter(result_df2)

    # Save the result to a CSV file
    result_path = os.path.join(resultsDir, f'{i}_{len(result_df2)}_analogs.csv')
    result_df2.to_csv(result_path, index=False)

    return result_df2, all_products

def reactBackTogether(analogs_reactant1: pd.DataFrame, analogs_reactant2: pd.DataFrame, ori_reaction, resultsDir):
    if ori_reaction is None:
        print('Must provide reaction to react back together.')
        # TODO: Could provide a way to react back together with all other reactions in the database.
        return None
    else:
        # Find all combinations of reactants that can be reacted back together with original reaction.
        # TODO: cycle through all reactions that you have SMARTS for. Find if reactant is True for reactant1 or reactant2.
        df_combined_results, all_products = findReactionCombinations(ori_reaction, analogs_reactant1, analogs_reactant2, resultsDir)

        # TODO: Perform structural scoring on the products to reference molecule.
        ori_mol = Chem.MolFromSmiles('CC1CCC(=CC(=O)NCCN2CCCCC2)C1')
        # structuralScores, similarMols = embedAndScoreSimilarSmiles(ori_mol, all_products)

    return df_combined_results

def searchReactantAnalogues(originalMol, reactant1, reactant2, ori_expansionAtomIdx=None,
                            similarityThr=config.SIMILARITY_SEARCH_THR,
                            structuralThr=config.STRUCTURAL_SCORE_THR, ori_reaction=None, resultsDir="./test/output2"):
    """
    Find analogues keeping the atoms involved in the reaction fixed.
    Author: Kate Fieseler

    :param originalMol: original molecule
    :param reactant1:
    :param reactant2:
    :param ori_expansionAtomIdx:
    :param similarityThr:
    :param structuralThr:
    :param ori_reaction: original reaction that links both reactants
    :param resultsDir:
    :return:
    """
    print('This is the original molecule:', Chem.MolToSmiles(originalMol))
    print('This is the first reactant:', Chem.MolToSmiles(reactant1))
    print('This is the second reactant:', Chem.MolToSmiles(reactant2))

    # IMPORTANT: ADDING HYDROGENS TO ALL MOLECULES
    Chem.AddHs(originalMol)
    Chem.AddHs(reactant1)
    Chem.AddHs(reactant2)

    # Find attachment points on reactants
    # Returns a list of tuples, (attachmentIdx_whole, attachmentIdx_subMol)
    # Find attachment points on reactants given reaction name
    reactant1_attachmentIdxs = find_attachmentIdxs_fromMolAndSubstruct(originalMol, reactant1, reaction=ori_reaction)
    reactant2_attachmentIdxs = find_attachmentIdxs_fromMolAndSubstruct(originalMol, reactant2, reaction=ori_reaction)

    print('reactant1_attachmentIdxs:', reactant1_attachmentIdxs)
    print('reactant2_attachmentIdxs:', reactant2_attachmentIdxs)

    if len(reactant1_attachmentIdxs) == 0 or len(reactant2_attachmentIdxs) == 0:
        print(f'No attachment points found for reactants... Stopping for this reaction {ori_reaction}')
        return None

    # Find atoms involved in the reaction
    # SMARTS pattern
    if ori_reaction is None:
        print('No reaction provided, need reaction name to find atoms involved in the reaction')
        return NotImplementedError  # TODO: Implement this!
    else:
        print('This is the original reaction:', ori_reaction)
        # Find atoms involved in the reaction
        # TODO: SEND THIS METHOD TO GRAVEYARD
        # react1_atoms = findReactionAtoms(originalMol, reactant1, ori_reaction, resultsDir, 1, reactant1_attachmentIdxs)
        # react2_atoms = findReactionAtoms(originalMol, reactant2, ori_reaction, resultsDir, 2, reactant2_attachmentIdxs)

        (r1_smarts, react1_atoms), (r2_smarts, react2_atoms) = findReactionAtoms_bothReactants(originalMol, reactant1, reactant2, ori_reaction, resultsDir, reactant1_attachmentIdxs, reactant2_attachmentIdxs)

    # Perform analogue search of reactants, keeping attachment point fixed
    analogs_react1 = processSubMolReactant(reactant1, ori_reaction, r1_smarts, react1_atoms, similarityThr, structuralThr, reactant_Num=1,
                          resultsDir=resultsDir)  # are the reacting atoms important here if indexing will change depending on the molecule?
    analogs_react2 = processSubMolReactant(reactant2, ori_reaction, r2_smarts, react2_atoms, similarityThr, structuralThr, reactant_Num=2,
                                           resultsDir=resultsDir)

    # Perform reaction search of reactant analogues (including PAINS filter)
    resultsdf = reactBackTogether(analogs_react1, analogs_react2, ori_reaction, resultsDir=resultsDir)

    # TODO: Do Fragmenstein placement within pipeline, adding failures to minimize to dataframe
    return resultsdf
