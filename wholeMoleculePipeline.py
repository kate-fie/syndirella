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
from chemUtils.distances import findCloseAtomIdxsFromRefAtomIdxs
from chemUtils.substructure import splitExtendedMolBySubstruct
from rdkit import Chem
from rdkit.Chem.PandasTools import WriteSDF
from rdkit.rdBase import BlockLogs

from config import config
from constants import fromReactionFullNameToReactantName, REACTIONS_NAMES
from databaseSearch.databaseSearch import perform_database_search
from embed import embedMolUsingRefs
from engine import myMap
from reactionSmarts import checkReactionSmartInAttachment, PATTERN_REACTANT1, PATTERN_REACTANT2, runReactionInAttachment
from scorer import scorePairOfMols


def getFunction_checkStructuralMatch(originalMol):
    def analyzeOneSmi(simMol):
        if not simMol:
            return None
        logsManager = BlockLogs()
        simEmbedded = embedMolUsingRefs(simMol, originalMol)
        del logsManager
        if simEmbedded is None:
            return None
        else:
            return scorePairOfMols(originalMol, simEmbedded)[0]

    return analyzeOneSmi

def embedAndScoreSimilarSmiles(refMol:Chem.Mol, smilesList) -> Tuple[Iterable[float], Iterable[Chem.Mol]]:
    """

    :param refMol: The reference mol to embed the smiles
    :param smilesList: A list of smiles
    :return: a list of tuples (structScore, smiles, similarity)
    """

    similar_mols = myMap(Chem.MolFromSmiles, smilesList)
    structuralScore = myMap(getFunction_checkStructuralMatch(refMol), similar_mols)
    return structuralScore, similar_mols

def matchMolToReactan(df, reactantNum):
    if reactantNum == 1:
        valid_patterns = PATTERN_REACTANT1
    elif reactantNum == 2:
        valid_patterns = PATTERN_REACTANT2
    else:
        raise ValueError("Only 1 ad 2 valid for reactant_num")
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

def createAnaloguesDf(smisDict, structuralScores, mols):

    similarities, metadatas = zip(* smisDict.values())
    result = OrderedDict(structuralScore=structuralScores, similarity=similarities, smi=smisDict.keys(),
                         metadata=metadatas, mol=mols)
    result = pd.DataFrame(result)
    result.dropna(axis=0, inplace=True)
    result.sort_values(by="structuralScore", inplace=True)
    result.reset_index(drop=True)
    return result

def processSubMol(subMol, attachmentIdx, similarityThr, structuralThr):

    smisDict = perform_database_search(subMol, thr=similarityThr)

    structuralScores, similarMols = embedAndScoreSimilarSmiles(subMol, smisDict.keys())

    df = createAnaloguesDf(smisDict, structuralScores, similarMols)

    equiv_attachmentIdxs = myMap(lambda mol: findCloseAtomIdxsFromRefAtomIdxs(mol, subMol, attachmentIdx), df["mol"])

    df.loc[:, "attachmentIdx"] = equiv_attachmentIdxs
    df = df[df["attachmentIdx"].notna()]
    if structuralThr is not None:
        df = df[df["structuralScore"] > structuralThr]
    df = matchMolToReactan(df, reactantNum=1)
    df = matchMolToReactan(df, reactantNum=2)
    df = df.reset_index(drop=True)
    return df

def searchDirectAnalogues(mol, similarityThr=config.SIMILARITY_SEARCH_THR, structuralThr=config.STRUCTURAL_SCORE_THR,
                          resultsDir="./output"):
    smisDict = perform_database_search(mol, thr=similarityThr)
    structuralScores, similarMols = embedAndScoreSimilarSmiles(mol,  smisDict.keys())
    results_df = createAnaloguesDf(smisDict, structuralScores, similarMols)
    #TODO: write results
    if structuralThr is not None:
        results_df = results_df[results_df["structuralScore"] > structuralThr]

    if resultsDir:
        sdfOutFname = os.path.join(resultsDir, "results.sdf")
        WriteSDF(results_df, sdfOutFname, molColName="mol", properties=results_df.columns[:-1])
        csvOutFname = os.path.join(resultsDir, "results.csv")
        results_df[[elem for elem in results_df.columns if elem != 'mols']].to_csv(csvOutFname, index=False)
    return results_df

def searchPicewiseAnalogues(originalMol, proposedMol, ori_expansionAtomIdx=None, similarityThr=config.SIMILARITY_SEARCH_THR,
                            structuralThr=config.STRUCTURAL_SCORE_THR, resultsDir="./output"):

    ori, grown, modified_attachmentIdx = splitExtendedMolBySubstruct(proposedMol, originalMol)
    (oriMol, ori_attachmentIdx) = ori
    (grownPart, grown_attachmentIdx) = grown

    # print(grownPart, grown_attachmentIdx)
    # from chemUtils.visualization import plotMols
    # plotMols([proposedMol, originalMol, oriMol, grownPart])

    if ori_expansionAtomIdx:
        assert ori_expansionAtomIdx == ori_attachmentIdx, "Error, proposedMol split didn't matched"

    ori_df = processSubMol(oriMol, ori_attachmentIdx, similarityThr, structuralThr)
    grown_df = processSubMol(grownPart, grown_attachmentIdx, similarityThr, structuralThr)

    if len(ori_df) == 0:
        error_msg = "Error, no valid analogues for original mol found!!"
        with open( os.path.join(resultsDir, "ERROR.txt"), "w") as f:
            f.write(error_msg)
        raise Exception(error_msg)

    if len(grown_df) == 0:
        error_msg = "Error, no valid analogues for attachment mol found!!"
        with open( os.path.join(resultsDir, "ERROR.txt"), "w") as f:
            f.write(error_msg)
        raise Exception(error_msg)

    if resultsDir:
        Chem.MolToMolFile(originalMol, os.path.join(resultsDir, "originalMol.mol"))
        Chem.MolToMolFile(proposedMol, os.path.join(resultsDir, "proposedMol.mol"))
        WriteSDF(ori_df, os.path.join(resultsDir, "ori_df.sdf"), molColName="mol", properties=ori_df.columns[:3])
        WriteSDF(grown_df, os.path.join(resultsDir, "grown_df.sdf"), molColName="mol", properties=grown_df.columns[:3])
        ori_df[[elem for elem in ori_df.columns if elem != 'mol']].to_csv( os.path.join(resultsDir, "ori_df.csv"), index=False)
        grown_df[[elem for elem in grown_df.columns if elem != 'mol']].to_csv( os.path.join(resultsDir, "grown_df.csv"), index=False)

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
                                                           reaction_products_similOri_similGrown_oriReact_grownReact)
        iter_score_simOri_simGrown_oriReact_grownReact = filter(None.__ne__, iter_score_simOri_simGrown_oriReact_grownReact)
        iter_score_simOri_simGrown_oriReact_grownReact = list(iter_score_simOri_simGrown_oriReact_grownReact)
        if len(iter_score_simOri_simGrown_oriReact_grownReact) == 0:
            #TODO: use logging
            print("No compounds found for %s"%reactionName)
            continue
        scores, products, similOri, similGrown, oriReact, grownReact  = zip(*iter_score_simOri_simGrown_oriReact_grownReact)

        oriSmiles = map(Chem.MolToSmiles, oriReact)
        grownSmiles = map(Chem.MolToSmiles, grownReact)

        results_df = pd.DataFrame(dict(structScore=scores, similarityFragment=similOri, similarityAttachment=similGrown,
                                       reactionName=reactionName, products=products, fragmentSMI=oriSmiles,
                                       attachmentSMI=grownSmiles))
        resultDfs_list.append(results_df)

    if len(resultDfs_list) == 0:
        error_msg = "Error, no valid products found!!"
        with open( os.path.join(resultsDir, "ERROR.txt"), "w") as f:
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
