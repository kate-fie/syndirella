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
from utils import input

def find_rows_to_process(df, batch_num, batch_size):
    # Get rows to be processed based on batch_num, batch_size or row argument
    if args.row is not None:
        rows_to_process = [args.row]
    elif args.batch is not None:
        start_idx, end_idx = calculate_batch_indices(len(df), args.batch_size, args.batch)
        rows_to_process = range(start_idx, end_idx)
    else:
        rows_to_process = df.index
    return rows_to_process

def calculate_batch_indices(total_rows, batch_size, batch_num):
    """
    Calculate the start and end indices for a batch.
    """
    start_idx = (batch_num - 1) * batch_size
    end_idx = min(start_idx + batch_size, total_rows)
    return start_idx, end_idx

def searchAnalogues(df, results_dir, superstructure, step_num, rows_to_process=1):
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
    # Convert NaN values to None
    df = df.applymap(lambda x: None if pd.isna(x) else x)

    for row_num in rows_to_process:
        row = df.iloc[row_num]
        if 'Done_Time' in df.columns:
            if row['Done_Time'] != None:
                continue
        if args.row is not None:
            if row_num != args.row:
                continue
        step = row['num_steps']
        """
        This part makes directories for each route and csv and runs the whole molecule pipeline
        
        INPUT TO WHOLE MOLECULE PIPELINE:
            ori_mol: RDKit mol object of the product
            reactant1_mol: RDKit mol object of the first reactant
            reactant2_mol: RDKit mol object of the second reactant
            ori_reaction: Reaction name to produce product
            resultsDirs: Path to save xtra results
            output_name: Name of the output csv (not full path)
        """
        ori_mol = Chem.MolFromSmiles((row['smiles']))
        if row['num_steps'] == 1:
            reaction_name = ast.literal_eval(row['rxn_order_first_to_last'])[0].replace(' ', '_')
        if reaction_name not in REACTIONS_NAMES:
            print(f"Do not have SMARTS for this reaction: {reaction_name}\n "
                  f"Please provide the SMARTS.\n"
                  f"Skipping {row['smiles']}...\n")
            continue
        else:
            reactants = ast.literal_eval(row['reactants'])[0] # CHANGE FOR OTHER STEPS
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
                xtra_reaction_dir_name = f"{results_dir}/{dir_name}/xtra_results_{step}/"
                os.makedirs(xtra_reaction_dir_name, exist_ok=True)
                xtra_reaction_dir_names.append(xtra_reaction_dir_name)
                output_name = f"{dir_name}_{step}_step"
            try:
                results = searchReactantAnalogues(ori_mol, reactant1_mol, reactant2_mol, step_num=step, ori_reaction=reaction_name,
                                              resultsDirs=xtra_reaction_dir_names, output_name=output_name, struct_score=False)
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
                results.to_csv(os.path.join(reaction_dir_names[i], f"{output_name}_{num_analogs}_1_of_2.csv"), index=False)
            print(results)
    if "OUTPUT" in args.input_csv:
        df.to_csv(args.input_csv, index=True)
    else:
        df.to_csv(f'{args.input_csv.split(".")[0]}_OUTPUT.csv', index=False)

    return os.path.join(reaction_dir_names[i], f"{output_name}_{num_analogs}.csv")

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
    parser.add_argument('-n', '--num_steps', help='Specify the number of steps in the route.', type=int, default=1)
    parser.add_argument('--test', help='Specify if you want to test the script.', action="store_true")
    parser.add_argument('--step1', help=".csv of the first step of multi step routes. Only needed for more than 1 step route.", type=str)


    args = parser.parse_args()

    if args.test and args.num_steps==2:
        # Load the smiles into dataframe
        df = pd.read_csv(args.input_csv, sep=',', header=0, index_col=0)
        # Get rows to be processed based on batch_num, batch_size or row argument
        rows_to_process = find_rows_to_process(df, args.batch, args.batch_size)
        searchAnalogues(df, args.results_dir, args.superstructure, step_num=2, rows_to_process=rows_to_process)

    if args.test:
        if args.step1 is not None:
            step1_df = pd.read_csv(args.step1, sep=',', header=0, index_col=0)
        else:
            df = pd.read_csv(args.input_csv, sep=',', header=0, index_col=0)
            step1_df = input.editstep1(df)
            step1_df.to_csv(f'{args.input_csv.split(".")[0]}_1st_step.csv', index=True)
        # df_edit = editstep1forstep2(df)
        searchAnalogues(df_edit, args.results_dir, args.superstructure, step_num=2)

    # Load the smiles into dataframe
    df = pd.read_csv(args.input_csv, sep=',', header=0, index_col=0)
    # Get rows to be processed based on batch_num, batch_size or row argument
    rows_to_process = find_rows_to_process(df, args.batch, args.batch_size)
    searchAnalogues(df, args.results_dir, args.superstructure, step_num=1, rows_to_process=rows_to_process)

    if args.num_steps == 2:
        df = pd.read_csv(args.input_csv, sep=',', header=0, index_col=0)
        df = input.editstep1(df)
        df.to_csv(f'{args.input_csv}_1st_step', index=True)
        rows_to_process = find_rows_to_process(df, args.batch, args.batch_size)
        step1_output = searchAnalogues(df, args.results_dir, args.superstructure, step_num=1)
        df = pd.read_csv(step1_output, sep=',', header=0, index_col=0)
        df_edit = editstep1forstep2(df)
        searchAnalogues(df_edit, args.results_dir, args.superstructure, step_num=2)



