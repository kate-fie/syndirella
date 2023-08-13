#!/usr/bin/env python3

"""
Created on Wed Oct 27 13:26:19 2021

Good lucking finding helpful analogues!

Author: Kate Fieseler
"""
import ast
import os
import sys
import argparse
import csv
import pandas as pd
from rdkit import Chem

sys.path.append('/Users/kate_fieseler/PycharmProjects/chemUtils')
import chemUtils

from config import config
from wholeMoleculePipeline import searchDirectAnalogues, searchPicewiseAnalogues, searchReactantAnalogues
from constants import REACTIONS_NAMES

def searchAnalogues(df, results_dir, superstructure):
    '''
    :param df: dataframe with columns SMILES, Reaction_name, Reactants
    :param results_dir: directory to save results
    :param superstructure: if True, perform superstructure search
    :return:
    '''
    # Create output directory
    try:
        os.mkdir(results_dir)
    except IOError:
        pass
    for i, row in df.iterrows():
        ori_mol = Chem.MolFromSmiles((row['SMILES']))
        reaction_name = row['Reaction_name'].replace(' ', '_')
        if reaction_name not in REACTIONS_NAMES:
            print(f"Do not have SMARTS for this reaction: {reaction_name}\n "
                  f"Please provide the SMARTS.\n"
                  f"Skipping {row['SMILES']}...\n")
            continue
        else:
            reactants = ast.literal_eval(row['Reactants'])
            print(reactants[0])
            reactant1_mol = Chem.MolFromSmiles(reactants[0])
            reactant2_mol = Chem.MolFromSmiles(reactants[1])
            print(reactant1_mol.GetNumAtoms())
            if reactant1_mol.GetNumAtoms() < 5 or reactant2_mol.GetNumAtoms() < 5:
                print("One reactant has very little atoms, the full pipeline search will not be performed since that "
                      "functionality is not implemented yet...\n")
                continue
            results = searchReactantAnalogues(ori_mol, reactant1_mol, reactant2_mol, ori_reaction=reaction_name, resultsDir=results_dir)
            results.to_csv(os.path.join(results_dir, f"{reaction_name}_{i}.csv"), index=False)
            print(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve synthesizable analogues from csv file of product SMILES, reactants, and reaction name.")
    parser.add_argument('-i', '--input_csv', type=str,
                        help=('Path to the input CSV file. The expected CSV structure is:\n'
                              'SMILES (str) - SMILES of product\n'
                              'Reaction_name (str) - Reaction name to produce product\n'
                              'Reactants (tuple) - Reactants listed in tuples\n'
                              '...\n'))
    parser.add_argument('-r', "--results_dir", help="Directory for the results", required=True)
    parser.add_argument('-u', "--superstructure", help='if performing a superstructure search', action="store_true")

    args = parser.parse_args()

    # TODO: Could parallelize search if searching for many SMILES
    # Load the smiles into dataframe
    df = pd.read_csv(args.input_csv)
    searchAnalogues(df, args.results_dir, args.superstructure)

