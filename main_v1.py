#!/usr/bin/env python3
# -*coding: utf-8 -*-
"""
Created on Wed Oct 27 13:59:55 2021

@author: ruben
"""
import os

from chemUtils.io import molFromFile, generatorFromFile

from config import config
from wholeMoleculePipeline import searchDirectAnalogues, searchPicewiseAnalogues


def cli_findAnalogues(initialMolFname: str, finalMolFname: str,
                      structuralScoreThr: float=config.STRUCTURAL_SCORE_THR,
                      similarityThr: float=config.SIMILARITY_SEARCH_THR,
                      outputDir: str = "./output", skipDirectAnalogues: bool = False,
                      skipBuildingbBlockAnalogues: bool = False,
                      ncpus: int = 1):
    '''

    :param initialMolFname: The original molecule (e.g. fragment)
    :param finalMolFname: The modified molecule derived from initialMol (e.g., version grown from initialMol)
    :param structuralScoreThr: The threshold for the structural score to select candidates (keep val>thr)
    :param similarityThr: The threshold for the similarity search (keep val>thr)
    :param outputDir: The directory where output will be generated
    :param skipDirectAnalogues: If True, do not performs regular similarity search
    :param skipBuildingbBlockAnalogues: If True, do not performs picewise similarity search
    :param ncpus: Number of CPUs for parallel execution
    :return:
    '''

    try:
        os.mkdir(outputDir)
    except IOError:
        pass

    initialMol = molFromFile(initialMolFname)
    finalMols = generatorFromFile(finalMolFname)

    for i, finalMol in enumerate(finalMols):
        if finalMol.HasProp("_Name"):
            molId = finalMol.GetProp("_Name")
            if molId =="":
                molId = str(i)
        else:
            molId = str(i)
        currentOutDir = os.path.join(outputDir, molId)
        try:
            os.mkdir(currentOutDir)
        except IOError:
            pass
        if not skipDirectAnalogues:
            outputSimilarityDir = os.path.join(currentOutDir, "similaritySearch")
            try:
                os.mkdir(outputSimilarityDir)
            except IOError:
                pass
            searchDirectAnalogues(finalMol, similarityThr=similarityThr, structuralThr=structuralScoreThr,
                                  resultsDir=outputSimilarityDir)

        if not skipBuildingbBlockAnalogues:
            outputSimilarBBsDir = os.path.join(currentOutDir, "similarBBlocks")
            try:
                os.mkdir(outputSimilarBBsDir)
            except IOError:
                pass
            searchPicewiseAnalogues(initialMol, finalMol, similarityThr=similarityThr,
                                    structuralThr=structuralScoreThr, resultsDir=outputSimilarBBsDir)

if __name__ == "__main__":

    from argParseFromDoc import get_parser_from_function

    parser = get_parser_from_function(cli_findAnalogues)
    args = parser.parse_args()
    print(cli_findAnalogues(**vars(args)))

'''
export PYTHONPATH=/home/sanchezg/sideProjects/chemUtil
python main.py --initialMolFname /home/sanchezg/oxford/myProjects/spend_65K_jun/round1_elaborating/dpp11/strife/input/bestCatalogue.sdf --finalMolFname /home/sanchezg/oxford/myProjects/spend_65K_jun/round1_elaborating/dpp11/strife/output/countsElabsMulti.sdf

'''