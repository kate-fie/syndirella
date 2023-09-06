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
from wholeMoleculePipeline import searchReactantAnalogues
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
        if args.row is not None:
            if i != args.row:
                continue
        ori_mol = Chem.MolFromSmiles((row['smiles']))
        if row['num_steps'] == 1:
            reaction_name = ast.literal_eval(row['rxn_order_first_to_last'].replace(' ', '_'))[0]
        if reaction_name not in REACTIONS_NAMES:
            print(f"Do not have SMARTS for this reaction: {reaction_name}\n "
                  f"Please provide the SMARTS.\n"
                  f"Skipping {row['smiles']}...\n")
            continue
        else:
            if row['num_steps'] == 1:
                reactants = ast.literal_eval(row['reactants'])[0] # for single step reactions
            print(reactants[0])
            reactant1_mol = Chem.MolFromSmiles(reactants[0])
            reactant2_mol = Chem.MolFromSmiles(reactants[1])
            if reactant1_mol.GetNumAtoms() < 5 or reactant2_mol.GetNumAtoms() < 5:
                print("One reactant has very little atoms, the full pipeline search will not be performed since that "
                      "functionality is not implemented yet...\n")
                continue
            dir_name = row['dir_name']
            reaction_dir_name = f"{results_dir}/{dir_name}/{reaction_name}_{Chem.MolToSmiles(ori_mol)}"
            os.makedirs(reaction_dir_name, exist_ok=True)
            # TODO: Add output name column to csv
            results = searchReactantAnalogues(ori_mol, reactant1_mol, reactant2_mol, ori_reaction=reaction_name,
                                              resultsDir=reaction_dir_name, output_name=dir_name+'-')
            if results is None:
                print("No results found for this molecule.\n")
                continue
            results.to_csv(os.path.join(results_dir, f"{reaction_name}_{i}.csv"), index=False)
            print(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve synthesizable analogues from csv file of product SMILES, reactants, and reaction name.")
    parser.add_argument('-i', '--input_csv', type=str,
                        help=('Path to the input CSV file. The expected CSV structure is:\n'
                              'SMILES (str) - SMILES of product\n'
                              'dir_name (str) - Name of the directory to save results. Usually target ID.\n'
                              'num_steps (int) - Number of steps in the route\n'
                              'rxn_order_first_to_last (list(str)) - Reaction name to produce product\n'
                              'reactants (list(tuple)) - Reactants listed in tuples\n'
                              '...\n'))
    parser.add_argument('-r', "--results_dir", help="Directory for the results", required=True)
    parser.add_argument('-u', "--superstructure", help='if performing a superstructure search', action="store_true")
    parser.add_argument('-b', '--row', help='specify row number to search. 0 is the first row below the header.', type=int)

    args = parser.parse_args()

    # TODO: Could parallelize search if searching for many SMILES
    # Load the smiles into dataframe
    df = pd.read_csv(args.input_csv)
    searchAnalogues(df, args.results_dir, args.superstructure)

