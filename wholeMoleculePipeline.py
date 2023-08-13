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
from chemUtils.geometry.molDistances import findCloseAtomIdxsFromRefAtomIdxs, find_close_atom_indices, find_substruct_atom_match
from chemUtils.substructure import splitExtendedMolBySubstruct
from chemUtils.substructure import find_attachmentIdxs_fromMolAndSubstruct
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.PandasTools import WriteSDF
from rdkit.rdBase import BlockLogs

from config import config
from constants import fromReactionFullNameToReactantName, REACTIONS_NAMES
from databaseSearch.databaseSearch import perform_database_search
from embed import embedMolUsingRefs
from engine import myMap
from reactionSmarts import checkReactionSmartInAttachment, PATTERN_REACTANT1, PATTERN_REACTANT2, \
    runReactionInAttachment, checkReactionSmartInAllAtoms, reaction_smarts
from scorer import scorePairOfMols

from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

import pickle
from itertools import chain

from PIL import Image, ImageDraw, ImageFont


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
        AllChem.EmbedMolecule(originalMol) # Lets try to embed the original mol, it worked
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


def matchMolToReactan(df, reactantNum):
    """
    Given a dataframe with a column "mol" and a column "attachmentIdx", this function checks if the molecule is a
    reactant of any the reactions in REACTIONS_NAMES.
    :param df:
    :param reactantNum:
    :return:
    """
    if reactantNum == 1:
        valid_patterns = PATTERN_REACTANT1
    elif reactantNum == 2:
        valid_patterns = PATTERN_REACTANT2  # Need to have two reactants, flip the order to find both combinations
    else:
        raise ValueError("Only 1 and 2 valid for reactant_num")
    matching_reactions = myMap(lambda mol_attachmentIdx:
                               checkReactionSmartInAttachment(*mol_attachmentIdx,
                                                              valid_patterns=valid_patterns),
                               zip(df["mol"], df["attachmentIdx"]))

    selected_reactions = {name: [] for name in REACTIONS_NAMES}
    for i, reactName_flag in enumerate(matching_reactions):
        for name, flag in reactName_flag:
            selected_reactions[name].append(flag)

    for name, vals in sorted(selected_reactions.items()):
        df.loc[:, fromReactionFullNameToReactantName(name, reactantNum)] = vals
    return df


def findReactionAtoms(ori_mol, reactant, reaction_name, results_dir):
    """
    Given a molecule and a reactant, this function finds the atoms that are involved in the reaction for the reactant.
    :param ori_mol:
    :param reactant:
    :param reaction_name:
    :return:
    """

    # TODO: Could provide option to add attachment index atoms when searching through the whole molecule with reaction smarts
    matched_reaction_atoms = checkReactionSmartInAllAtoms(reactant, reaction_name)
    if matched_reaction_atoms is None:
        print('No reaction atoms found for reactant', Chem.MolToSmiles(reactant), 'in reaction', reaction_name)
        return None
    matched_reaction_atoms = tuple(x for sub in matched_reaction_atoms for x in sub)  # This flattens the tuple
    save_molecule_with_indices(reactant, os.path.join(results_dir,
                                                      f"{reaction_name}_reactant{Chem.MolToSmiles(reactant)}.png"), matched_reaction_atoms)

    return matched_reaction_atoms


def createAnaloguesDf(smisDict, structuralScores, mols, mw_diff, mw, similarity):
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
    result.sort_values(by="structuralScore", inplace=True)
    result.reset_index(drop=True, inplace=True)
    return result


def processSubMol(subMol, attachmentIdx, similarityThr, structuralThr, resultsDir="./test/output"):
    """
    Given a substructure and an attachment point, this function searches for analogues of the substructure

    :param subMol:
    :param attachmentIdx:
    :param similarityThr:
    :param structuralThr:
    :return:
    """

    print('This is the substructure:', Chem.MolToSmiles(subMol))
    print('This is the attachment point:', attachmentIdx)

    smisDict = perform_database_search(subMol, thr=similarityThr)

    # Save smisDict to a file, to show structures easily
    out_dict = os.path.join(resultsDir, 'smisDict.pickle')
    with open(out_dict, 'wb') as f:
        pickle.dump(smisDict, f)

    structuralScores, similarMols = embedAndScoreSimilarSmiles(subMol, smisDict.keys())

    df = createAnaloguesDf(smisDict, structuralScores, similarMols)

    equiv_attachmentIdxs = myMap(lambda mol: findCloseAtomIdxsFromRefAtomIdxs(mol, subMol, attachmentIdx), df[
        "mol"])  # TODO: understand this, I think it finds the same attachment point as original molecule

    df.loc[:, "attachmentIdx"] = equiv_attachmentIdxs
    df = df[df["attachmentIdx"].notna()]
    if structuralThr is not None:
        df = df[df["structuralScore"] > structuralThr]
    df = matchMolToReactan(df, reactantNum=1)
    df = matchMolToReactan(df, reactantNum=2)
    df = df.reset_index(drop=True)
    return df


def processSubMolReactant(reactant, reaction_atoms, similarityThr, structuralThr, resultsDir, substructure=True):
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

    smisDict = perform_database_search(reactant, substructure=substructure, thr=similarityThr) # substructure is the search for superstructure

    # Save smisDict to a file, to show structures easily
    out_dict = os.path.join(resultsDir, 'smisDict.pickle')
    with open(out_dict, 'wb') as f:
        pickle.dump(smisDict, f)

    structuralScores, similarMols, mw_diff, mw = embedAndScoreSimilarSmiles(reactant, smisDict.keys())

    similarity = True
    if substructure is True: similarity = False

    df = createAnaloguesDf(smisDict, structuralScores, similarMols, mw_diff, mw, similarity)
    df.to_csv(os.path.join(resultsDir, f'analogues_preprocess_{Chem.MolToSmiles(reactant)}.csv'))

    # Search through elaborated compounds for the reaction atoms as substructure search
    equiv_attachmentIdxs = myMap(lambda mol: find_substruct_atom_match(mol, reactant, reaction_atoms), df[
        "mol"])

    df.loc[:, "attachmentIdx"] = equiv_attachmentIdxs
    df = df[df["attachmentIdx"].notna()]
    print('This is the length of the dataframe after substructure search:', len(df))
    if structuralThr is not None:
        df = df[df["structuralScore"] > structuralThr]
    df = matchMolToReactan(df, reactantNum=1)
    df = matchMolToReactan(df, reactantNum=2)
    df = df.reset_index(drop=True)
    df.to_csv(os.path.join(resultsDir, f'analogues_substruct_{Chem.MolToSmiles(reactant)}.csv'))
    return df

import pandas as pd

import pandas as pd
from itertools import product

def apply_reaction(row, reaction_name):
    """
    Applies a reaction to a row of a dataframe. Need to make sure that the reaction SMARTS is flipped so both sides of
    the reaction are tried out. But need to make sure that the reactants are in the right order like how they were originally discovered.
    Since functional groups decorating a ring do not want to be reacted together even if they match the smarts pattern.

    :param row:
    :param reaction_name:
    :return: product SMILES
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

    filtered_df1 = filtered_df1[['structuralScore', 'mw', 'smi', 'metadata', 'attachmentIdx']]
    filtered_df2 = filtered_df2[['structuralScore', 'mw', 'smi', 'metadata', 'attachmentIdx']]

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
        'attachmentIdx_reactant1': 'reactant1_attachmentIdx',
        'smi_reactant2': 'reactant2_smi',
        'structuralScore_reactant2': 'reactant2_structuralScore',
        'mw_reactant2': 'reactant2_mw',
        'metadata_reactant2': 'reactant2_metadata',
        'attachmentIdx_reactant2': 'reactant2_attachmentIdx'
    }, inplace=True)

    # Add a column of the reacted together SMILES
    result_df['product_smi'] = result_df.apply(apply_reaction, reaction_name=i, axis=1)

    # Extracting the product_smi column and flattening the lists
    all_products = list(set(chain.from_iterable(result_df['product_smi'])))

    # Save the result to a CSV file
    result_path = os.path.join(resultsDir, f'{i}_{len(result_df)}_analogs.csv')
    result_df.to_csv(result_path, index=False)

    return result_df, all_products

def reactBackTogether(analogs_reactant1, analogs_reactant2, ori_reaction, resultsDir):
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

def searchDirectAnalogues(mol, similarityThr=config.SIMILARITY_SEARCH_THR, structuralThr=config.STRUCTURAL_SCORE_THR,
                          resultsDir="./output"):
    smisDict = perform_database_search(mol, thr=similarityThr)
    structuralScores, similarMols = embedAndScoreSimilarSmiles(mol, smisDict.keys())
    results_df = createAnaloguesDf(smisDict, structuralScores, similarMols)
    # TODO: write results
    if structuralThr is not None:
        results_df = results_df[results_df["structuralScore"] > structuralThr]

    if resultsDir:
        sdfOutFname = os.path.join(resultsDir, "results.sdf")
        WriteSDF(results_df, sdfOutFname, molColName="mol", properties=results_df.columns[:-1])
        csvOutFname = os.path.join(resultsDir, "results.csv")
        results_df[[elem for elem in results_df.columns if elem != 'mols']].to_csv(csvOutFname, index=False)
    return results_df


def searchPicewiseAnalogues(originalMol, proposedMol, ori_expansionAtomIdx=None,
                            similarityThr=config.SIMILARITY_SEARCH_THR,
                            structuralThr=config.STRUCTURAL_SCORE_THR, resultsDir="./test/output"):
    """
    Given an original mol and proposed mol (which is an extension of the original mol), this function splits the
    proposed mol in two parts: the original mol and the grown part. Then, it searches for analogues of the grown part
    and for analogues of the original mol. Finally, it combines the results of the two searches and returns a dataframe
    with the results.

    :param originalMol:
    :param proposedMol:
    :param ori_expansionAtomIdx:
    :param similarityThr:
    :param structuralThr:
    :param resultsDir:
    :return:
    """

    # ori = original mol, grown = grown part of the proposed mol
    ori, grown, modified_attachmentIdx = splitExtendedMolBySubstruct(proposedMol, originalMol)
    (oriMol, ori_attachmentIdx) = ori
    (grownPart, grown_attachmentIdx) = grown

    # print(grownPart, grown_attachmentIdx)
    # from chemUtils.visualization import plotMols
    # plotMols([proposedMol, originalMol, oriMol, grownPart])

    if ori_expansionAtomIdx:
        assert ori_expansionAtomIdx == ori_attachmentIdx, "Error, proposedMol split didn't matched"

    print('This is the original mol:', Chem.MolToSmiles(oriMol))
    print('This is the grown part not involved in linking to original:', Chem.MolToSmiles(grownPart))

    # ori_df = analogues of the original mol, grown_df = analogues of the grown part
    ori_df = processSubMol(oriMol, ori_attachmentIdx, similarityThr, structuralThr)
    grown_df = processSubMol(grownPart, grown_attachmentIdx, similarityThr, structuralThr)

    if len(ori_df) == 0:
        error_msg = "Error, no valid analogues for original mol found!!"  # Getting this error
        with open(os.path.join(resultsDir, "ERROR.txt"), "w") as f:
            f.write(error_msg)
        raise Exception(error_msg)

    if len(grown_df) == 0:
        error_msg = "Error, no valid analogues for attachment mol found!!"
        with open(os.path.join(resultsDir, "ERROR.txt"), "w") as f:
            f.write(error_msg)
        raise Exception(error_msg)

    if resultsDir:
        Chem.MolToMolFile(originalMol, os.path.join(resultsDir, "originalMol.mol"))
        Chem.MolToMolFile(proposedMol, os.path.join(resultsDir, "proposedMol.mol"))
        WriteSDF(ori_df, os.path.join(resultsDir, "ori_df.sdf"), molColName="mol", properties=ori_df.columns[:3])
        WriteSDF(grown_df, os.path.join(resultsDir, "grown_df.sdf"), molColName="mol", properties=grown_df.columns[:3])
        ori_df[[elem for elem in ori_df.columns if elem != 'mol']].to_csv(os.path.join(resultsDir, "ori_df.csv"),
                                                                          index=False)
        grown_df[[elem for elem in grown_df.columns if elem != 'mol']].to_csv(os.path.join(resultsDir, "grown_df.csv"),
                                                                              index=False)

    resultDfs_list = []
    for reactionName in REACTIONS_NAMES:
        react1_name = fromReactionFullNameToReactantName(reactionName, 1)
        react2_name = fromReactionFullNameToReactantName(reactionName, 2)

        # TODO: check corner case in which no matching idx is found in np.where
        valid_ori_1 = np.where(ori_df[react1_name])[0]
        valid_grown_2 = np.where(grown_df[react2_name])[0]
        valid_ori_2 = np.where(ori_df[react2_name])[0]
        valid_grown_1 = np.where(grown_df[react1_name])[0]
        candidateIdxs = itertools.chain.from_iterable([
            itertools.product(valid_ori_1, valid_grown_2),
            itertools.product(valid_ori_2, valid_grown_1)])
        candidateIdxs = list(candidateIdxs)
        if len(candidateIdxs) == 0:
            continue

        def runReactants(r1Idx_r2Idx):

            oriIdx, grownIdx = r1Idx_r2Idx
            ori_reactant, ori_attachmentIdxs, ori_simil = ori_df.iloc[oriIdx, :][["mol", "attachmentIdx", "similarity"]]
            grown_reactant, grown_attachmentIdxs, grown_simil = grown_df.iloc[grownIdx, :][
                ["mol", "attachmentIdx", "similarity"]]

            candidate_products = runReactionInAttachment(ori_reactant, grown_reactant, reactionName, refMol=proposedMol,
                                                         ref_attachmentIdxs=modified_attachmentIdx)

            # from molPlot import plotMols
            # print(Chem.MolToSmiles(ori_reactant), Chem.MolToSmiles(grown_reactant), Chem.MolToSmiles(proposedMol))
            # print(reactionName, ori_attachmentIdxs, grown_attachmentIdxs, modified_attachmentIdx)
            # plotMols([ori_reactant, grown_reactant, proposedMol])

            return candidate_products, ori_simil, grown_simil, ori_reactant, grown_reactant

        reaction_products_similOri_similGrown_oriReact_grownReact = myMap(runReactants, candidateIdxs)

        def processProductList(productList_andMetadata):
            '''

            :param productList_andMetadata: A list with the first element being a productMol and the other elemnts metadata
            :return:
            '''

            listOfProducts = productList_andMetadata[0]
            if not listOfProducts:
                return None
            score_prod = []
            for product in listOfProducts:
                score = scorePairOfMols(proposedMol, product)[0]
                score_prod.append((score, product))
            score_prod = sorted(score_prod)
            return score_prod[-1] + productList_andMetadata[1:]

        iter_score_simOri_simGrown_oriReact_grownReact = myMap(processProductList,
                                                               reaction_products_similOri_similGrown_oriReact_grownReact)  # TODO: Not even getting to this line
        iter_score_simOri_simGrown_oriReact_grownReact = filter(None.__ne__,
                                                                iter_score_simOri_simGrown_oriReact_grownReact)
        iter_score_simOri_simGrown_oriReact_grownReact = list(iter_score_simOri_simGrown_oriReact_grownReact)
        if len(iter_score_simOri_simGrown_oriReact_grownReact) == 0:
            # TODO: use logging
            print("No compounds found for %s" % reactionName)
            continue
        scores, products, similOri, similGrown, oriReact, grownReact = zip(
            *iter_score_simOri_simGrown_oriReact_grownReact)

        oriSmiles = map(Chem.MolToSmiles, oriReact)
        grownSmiles = map(Chem.MolToSmiles, grownReact)

        results_df = pd.DataFrame(dict(structScore=scores, similarityFragment=similOri, similarityAttachment=similGrown,
                                       reactionName=reactionName, products=products, fragmentSMI=oriSmiles,
                                       attachmentSMI=grownSmiles))
        resultDfs_list.append(results_df)

    if len(resultDfs_list) == 0:
        error_msg = "Error, no valid products found!!"  # TODO: Figure out why getting this error...
        with open(os.path.join(resultsDir, "ERROR.txt"), "w") as f:
            f.write(error_msg)
        raise Exception(error_msg)

    results_df = pd.concat(resultDfs_list)
    results_df.loc[:, "productSMI"] = results_df.loc[:, "products"].map(Chem.MolToSmiles)
    results_df = results_df.dropna(axis=0).sort_values(by="structScore", ascending=False).reset_index(drop=True)

    if resultsDir:
        sdfOutFname = os.path.join(resultsDir, "results.sdf")
        WriteSDF(results_df, sdfOutFname, molColName="products", properties=results_df.columns[:3])
        csvOutFname = os.path.join(resultsDir, "results.csv")
        results_df[[elem for elem in results_df.columns if elem != 'products']].to_csv(csvOutFname, index=False)
    return results_df


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

    # Find attachment points on reactants
    reactant1_attachmentIdxs = find_attachmentIdxs_fromMolAndSubstruct(originalMol, reactant1)
    reactant2_attachmentIdxs = find_attachmentIdxs_fromMolAndSubstruct(originalMol, reactant2)

    print('reactant1_attachmentIdxs:', reactant1_attachmentIdxs)
    print('reactant2_attachmentIdxs:', reactant2_attachmentIdxs)

    # Find atoms involved in the reaction
    # SMARTS pattern
    if ori_reaction is None:
        print('No reaction provided, need reaction name to find atoms involved in the reaction')
        return NotImplementedError  # TODO: Implement this!
    else:
        print('This is the original reaction:', ori_reaction)
        # Find atoms involved in the reaction
        react1_atoms = findReactionAtoms(originalMol, reactant1, ori_reaction, resultsDir)
        react2_atoms = findReactionAtoms(originalMol, reactant2, ori_reaction, resultsDir)

    # Perform analogue search of reactants, keeping attachment point fixed
    analogs_react1 = processSubMolReactant(reactant2, react2_atoms, similarityThr, structuralThr,
                          resultsDir=resultsDir)  # are the reacting atoms important here if indexing will change depending on the molecule?
    analogs_react2 = processSubMolReactant(reactant1, react1_atoms, similarityThr, structuralThr, resultsDir=resultsDir)

    # Perform reaction search of reactant analogues
    resultsdf = reactBackTogether(analogs_react1, analogs_react2, ori_reaction, resultsDir=resultsDir)

    return resultsdf # TODO: Implement this!
