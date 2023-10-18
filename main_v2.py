#!/usr/bin/env python3

"""
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
import datetime
import numpy as np

sys.path.append('/Users/kate_fieseler/PycharmProjects/chemUtils')
import chemUtils

from config import config
from wholeMoleculePipeline import searchReactantAnalogues
from constants import REACTIONS_NAMES

def calculate_batch_indices(total_rows, batch_size, batch_num):
    """
    Calculate the start and end indices for a batch.
    """
    start_idx = (batch_num - 1) * batch_size
    end_idx = min(start_idx + batch_size, total_rows)
    return start_idx, end_idx

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
    if 'Done_Time' not in df.columns:
        df['Done_Time'] = np.nan
    if 'Analogs' not in df.columns:
        df['Analogs'] = np.nan
    # Get rows to be processed based on batch_num, batch_size or row argument
    if args.row is not None:
        rows_to_process = [args.row]
    elif args.batch is not None:
        start_idx, end_idx = calculate_batch_indices(len(df), args.batch_size, args.batch)
        rows_to_process = range(start_idx, end_idx)
    else:
        rows_to_process = df.index

    for row_num in rows_to_process:
        row = df.iloc[row_num]
        if args.row is not None:
            if row_num != args.row:
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
            dir_names = ast.literal_eval(row['dir_name'])
            # Remove colon from dir names
            dir_names = [dir_name.replace(':', '_') for dir_name in dir_names]
            reaction_dir_names = []
            xtra_reaction_dir_names = []
            for i in range(len(dir_names)):
                dir_name = dir_names[i]
                reaction_dir_name = f"{results_dir}/{dir_name}/"
                os.makedirs(reaction_dir_name, exist_ok=True)
                reaction_dir_names.append(reaction_dir_name)
                xtra_reaction_dir_name = f"{results_dir}/{dir_name}/xtra_results/"
                os.makedirs(xtra_reaction_dir_name, exist_ok=True)
                xtra_reaction_dir_names.append(xtra_reaction_dir_name)
                output_name = f"{dir_name}_{row['num_steps']}_step"

            try:
                results = searchReactantAnalogues(ori_mol, reactant1_mol, reactant2_mol, ori_reaction=reaction_name,
                                              resultsDirs=xtra_reaction_dir_names, output_name=output_name)
                # Set the Done_Time and Analogs columns for the current row
                df.at[row_num, 'Done_Time'] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                df.at[row_num, 'Analogs'] = len(results) if results is not None else 0
            except Exception as e:
                print(f"Error processing row {i}: {e}")
                continue

            if results is None:
                print("No results found for this molecule.\n")
                continue
            for i in range(len(reaction_dir_names)):
                num_analogs = len(results)
                results.to_csv(os.path.join(reaction_dir_names[i], f"{output_name}_{num_analogs}.csv"), index=False)
            print(results)
    if "OUTPUT" in args.input_csv:
        df.to_csv(args.input_csv, index=True)
    else:
        df.to_csv(f'{args.input_csv.split(".")[0]}_OUTPUT.csv', index=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve synthesizable analogues from csv file of product SMILES, reactants, and reaction name.")
    parser.add_argument('-i', '--input_csv', type=str,
                        help=('Path to the input CSV file. The expected CSV structure is:\n'
                              'SMILES (str) - SMILES of product\n'
                              'dir_name (list) - Name of the directories to save results. Usually target ID.\n'
                              'num_steps (int) - Number of steps in the route\n'
                              'rxn_order_first_to_last (list(str)) - Reaction name to produce product\n'
                              'reactants (list(tuple)) - Reactants listed in tuples\n'
                              '...\n'))
    parser.add_argument('-r', "--results_dir", help="Directory for the results", required=True)
    parser.add_argument('-u', "--superstructure", help='if performing a superstructure search', action="store_true")
    parser.add_argument('-b', '--row', help='specify row number to search. 0 is the first row below the header.', type=int)
    parser.add_argument('--batch', help='specify batch number to search. 1 is the first batch.', type=int)
    parser.add_argument('--batch_size', help='specify batch size. Default is 10.', type=int, default=10)
    parser.add_argument('--go', help='Specify if you want to keep going after starting a batch.', action="store_true")


    args = parser.parse_args()

    # TODO: Could parallelize search if searching for many SMILES
    # Load the smiles into dataframe
    df = pd.read_csv(args.input_csv, sep=',', header=0, index_col=0)
    searchAnalogues(df, args.results_dir, args.superstructure)

