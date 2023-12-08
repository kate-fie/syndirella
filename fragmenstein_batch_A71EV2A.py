#!/usr/bin/env python
"""
Author: Kate Fieseler
Created: 17 October 2023
Description: Batch script for running Fragmenstein
"""
import os
import sys
import argparse
import subprocess
import re

import pandas as pd
import datetime
import shutil
import traceback

from typing import List, Dict, Any, Optional

def config_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', required=True, help='Home directory where csvs of elabs live.')
    parser.add_argument('-t', required=True, help='Directory of templates to use. Name must contain fragalysis ID.')
    parser.add_argument('-i', required=True, help='SDF of all fragment hits for target.')
    parser.add_argument('--step', required=True, help='Step identifier to place elabs for. Must be present in csv name.')
    parser.add_argument('--batch_range', required=True,
                        help='Num atom difference range for placing. E.g., "0-3", "4-7", "8-11", "12+".')
    parser.add_argument('--exact_template', required=False, help='Template pdb to use for all elabs.')
    parser.add_argument('-o', required=False, help='Output csv to save logging results.')
    parser.add_argument('-p', required=False, help='Prefix for fragment names in sdf.')
    parser.add_argument('--n_cores', required=False, default=1, help='Number of cores to use.')
    parser.add_argument('-l', '--log_path', required=False, help='Path to the log directory.')
    parser.add_argument('--cutoff', action='store_true', required=False, help='Cutoff of 5,000 for placing elabs.')
    parser.add_argument('--wictor', action='store_true', required=False, help='Use if running on with Wictor in Fragmenstein')
    return parser


def extract_info(directory):
    """Extracts the cmpd_catalog, frag1, and frag2 from the directory name."""
    # Define a regular expression pattern to capture the cmpd_catalog, frag1, and frag2
    pattern = r'(.+?)_(x\d+_\d+[AB])_(x\d+_\d+[AB])'
    # Use re.search to find the first match in the input string
    match = re.search(pattern, directory)
    if match:
        # Extract the cmpd_catalog, frag1, and frag2 from the matched groups
        cmpd_catalog, frag1, frag2 = match.groups()
        return cmpd_catalog, frag1, frag2
    else:
        # Return None if the pattern is not found in the input string
        return None, None, None

def shorten_elabs_csv(elabs_csv, len):
    """Shortens the elabs csv to 10,000 rows."""
    df = pd.read_csv(elabs_csv, encoding='ISO-8859-1')
    df = df[:10000]
    new_path = elabs_csv.replace(f'.csv','') + "_10000.csv"
    df.to_csv(new_path, index=True)
    return new_path

def extract_molecule_from_sdf(sdf_content, molecule_name):
    """Extracts the molecule with the given name from the SDF content."""
    molecules = sdf_content.split("$$$$\n")
    for molecule in molecules:
        if molecule_name in molecule:
            return molecule + "$$$$\n"
    return None

def filter_sdf_by_names(sdf_content, molecule_names: list, sdf_prefix=None):
    """Filters the SDF content to retain only the molecules with the given names."""
    if sdf_prefix:
        molecule_names = [f"{sdf_prefix}{molecule_name}" for molecule_name in molecule_names]
    return "".join(extract_molecule_from_sdf(sdf_content, name) for name in molecule_names)

def find_template_pdb(elabs_csv, template_dir, kwargs):
    """Finds the template pdb as a column in the csv."""
    if kwargs['exact_template']:
        ref_pdb = kwargs['exact_template']
    else:
        df = pd.read_csv(elabs_csv, encoding='ISO-8859-1')
        if 'ref_pdb' in df.columns:
            ref_pdb = df['ref_pdb'][0]
        else:
            print("ref_pdb column not found in csv.")
            return None
    for root, dirs, files in os.walk(template_dir):
        for file in files:
            if ref_pdb in file and not file.startswith('.'): # ignore hidden
                return os.path.join(root, file)
    return None

def find_frags_sdf(sdf_content, root, directory, cmpd_catalog, frag1, frag2, sdf_prefix=None):
    """Finds the sdf of the frags for the given cmpd_catalog, frag1, and frag2."""
    frags_sdf = os.path.join(root, directory, f"{cmpd_catalog}_{frag1}_{frag2}_frags.sdf")
    if os.path.exists(frags_sdf) is False:
        frag_sdf_content = filter_sdf_by_names(sdf_content, [frag1, frag2], sdf_prefix=sdf_prefix)
        with open(frags_sdf, 'w') as f:
            f.write(frag_sdf_content)
    return frags_sdf

def find_elabs_csv(root, directory, step_identifier):
    """Finds the elabs csv for the given directory. Adds column 'hit_names' if needed."""
    for filename in os.listdir(os.path.join(root, directory)):
        if filename.endswith('.csv') and step_identifier in filename and not filename.startswith('.') \
                and directory in filename and 'batch' not in filename:
            csv_path = os.path.join(root, directory, filename)
            df = pd.read_csv(csv_path, encoding='ISO-8859-1')
            return csv_path, len(df)
    return None, None

# Function to search within lists in the DataFrame column
def search_in_list(row, dir_name):
    return dir_name in row

def set_acceptable_value(row, matches_index, num_accept_elabs):
    if not row.empty and row.name in matches_index:
        row['placed_acceptable'] = num_accept_elabs  # assuming you have this from the output csv processing

def update_placed_acceptable(df, directory, num_accept_value):
    print(f"APPENDING {directory} to all elabs csv.")
    # Ensure 'placed_acceptable' column exists
    if 'placed_acceptable' not in df.columns:
        df['placed_acceptable'] = None  # or pd.NA for pandas' native missing value type
    # Find the index of rows with the matching 'dir_name'
    matches_indices = df[df['dir_name'].apply(lambda x: directory in x)].index
    # Update the 'placed_acceptable' column for matched rows
    df.loc[matches_indices, 'placed_acceptable'] = num_accept_value
    return df

def find_p1_frag(frag1, frag2):
    """Finds the frag in the P1 site out of the 2."""
    frag1 = frag1.split("_")[0]
    frag2 = frag2.split("_")[0]
    if frag1 == frag2:
        return frag1
    else:
        return None

def run_batch(**kwargs):
    """
    Runs:
    fragmenstein laboratory place --input hits.sdf --template template.pdb --in-table test.csv

    Where:
    "--input" = fragments .sdf
    "--template" = apo desolv of first fragment in name
    "--in-table" = csv of elabs to place

    INPUT:
        kwargs['h']:home directory of all elabs organized in respective folders with name:
                    base_catalog_name_frag1_frag2 (ex. ENAMINE_PV-000782638925_x1052_1A_x1083_0A)
        kwargs['t']:directory of apo desolv templates to use. Each must contain fragID
                    (ex. D68EV3CPROA-x0102_1B_apo-desolv.pdb)
        kwargs['i']:sdf of ALL fragment hits.
        kwargs['o']:name with path to save summarizing csv to.
    """
    if not os.path.exists(kwargs['d']):
        print("Home directory does not exist.")
    if not os.path.exists(kwargs['t']):
        print("Template directory does not exist.")
    if not os.path.exists(kwargs['i']):
        print("SDF of all fragment hits does not exist.")

    sdf_file_path = kwargs['i']
    n_cores = kwargs['n_cores']
    with open(sdf_file_path, 'r') as f:
        sdf_content = f.read()

    batch_range = kwargs['batch_range']
    min_value, max_value = None, None

    # Parse the batch range
    if batch_range == "12+":
        min_value = 12
    else:
        min_value, max_value = map(int, batch_range.split('-'))

    for root, dirs, files in os.walk(kwargs['d']):
        for directory in dirs:
            print('DIRECTORY', directory)
            if directory == 'output':
                exit()
            if "xtra_results" in directory or "logs" in directory:
                exit()
            elabs_csv, len = find_elabs_csv(root, directory, step_identifier=kwargs['step'])
            if elabs_csv is None:
                print(f"Elabs csv not found for {directory}.")
                continue
            print(elabs_csv)
            frags_sdf = sdf_file_path
            print(frags_sdf)
            template_pdb = find_template_pdb(elabs_csv, kwargs['t'], kwargs)
            print(template_pdb)
            # Filter the DataFrame based on the batch range
            df = pd.read_csv(elabs_csv, encoding='ISO-8859-1')
            if max_value is None:  # For "12+" case
                batch_df = df[df['num_atom_difference'] >= min_value]
            else:
                batch_df = df[df['num_atom_difference'].between(min_value, max_value)]
            batch_csv_path = None
            if not batch_df.empty:
                batch_len = batch_df.shape[0]
                batch_csv_path = elabs_csv.replace('.csv', f'_batch_{batch_range}.csv')
                output_batch_csv_path = batch_csv_path.replace('.csv', f'_output.csv')
                batch_df.to_csv(batch_csv_path, index=False)
            if frags_sdf is None:
                print(f"Frags sdf not found for {directory}.")
            if template_pdb is None:
                print(f"Template pdb not found for {directory}.")
            if frags_sdf and template_pdb and batch_csv_path:
                os.chdir(os.path.join(root, directory)) # change to directory
                try:
                    if kwargs['cutoff'] and batch_len > 10000:
                        print(f"CUTTING {directory} because it has more than 10,000 elabs.")
                        batch_csv_path = shorten_elabs_csv(batch_csv_path, len)
                    if kwargs['wictor']:
                        print(f"PLACING {directory} WITH WICTOR.")
                        command = ["fragmenstein", "laboratory", "place", "--input", frags_sdf, "--template",
                                   template_pdb,
                                   "--in-table", batch_csv_path, "--output", output_batch_csv_path, "--cores", str(n_cores),
                                   "--verbose", "--victor", "Wictor"]
                        command_str = ' '.join(command)
                        print(f"Executing command: {command_str}")
                        subprocess.run(command)
                    else:
                        print(f"PLACING {directory}.")
                        command = ["fragmenstein", "laboratory", "place", "--input", frags_sdf, "--template", template_pdb,
                             "--in-table", batch_csv_path, "--output", output_batch_csv_path, "--cores", str(n_cores), "--verbose"]
                        command_str = ' '.join(command)
                        print(f"Executing command: {command_str}")
                        subprocess.run(command)
                    print("DELETING EXTRA FILES.")
                    output_dir = os.path.join(root, directory, 'output')
                    if not os.path.exists(output_dir):
                        print(f"Output directory does not exist for {directory}.")
                        continue
                    os.chdir(output_dir)  # change to output directory
                    suffixes_to_keep = ['.minimised.json', '.minimised.mol', '.csv']
                    for root1, dirs, files in os.walk(output_dir):
                        for directory2 in dirs:
                            directory2_path = os.path.join(output_dir, directory2)
                            if not os.path.exists(directory2_path):
                                print('THIS PATH DOES NOT EXIST', directory2_path)
                                continue
                            if 'base' in directory2:
                                print(f"Skipping {directory2} because it is the base compound.")
                                continue
                            for file in os.listdir(directory2_path):
                                if not any(file.endswith(suffix) for suffix in suffixes_to_keep):
                                    file_path = os.path.join(output_dir, directory2, file)
                                    try:
                                        if os.path.isfile(file_path) or os.path.islink(file_path):
                                            os.unlink(file_path)
                                            print(f"Deleted file: {file_path}")
                                        elif os.path.isdir(file_path):
                                            shutil.rmtree(file_path)
                                            print(f"Deleted directory: {file_path}")
                                    except Exception as e:
                                        print('Failed to delete %s. Reason: %s' % (file_path, e))
                except Exception as e:
                    tb = traceback.format_exc()
                    print(f"Error placing elabs for {directory}.")
                    print(tb)
                    continue
            else:
                continue

def main():
    parser = config_parser()
    # load
    settings: Dict[str, Any] = vars(parser.parse_args())
    # Redirect stdout and stderr to log file
    # Generate a timestamp string and append to the log file path
    if settings['log_path'] is None:
        # run
        try:
            print('RUNNING')
            run_batch(**settings)
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)
    else:
        # I don't think this logging functionality works
        if os.path.exists(settings['log_path']) is False:
            os.makedirs(settings['log_path'])
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        log_path_with_timestamp = f"{settings['log_path']}frag_place_{timestamp}.log"
        with open(log_path_with_timestamp, 'a') as f:
            sys.stdout = f
            sys.stderr = f
            print("Arguments passed to the script:")
            for key, value in settings.items():
                print(f"{key}: {value}")
            # run
            try:
                run_batch(**settings)
            except Exception as e:
                print(f"Error: {e}", file=sys.stderr)
                sys.exit(1)


if __name__ == "__main__":
    main()