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

import pandas as pd
import datetime
import shutil

from typing import List, Dict, Any, Optional

def config_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', required=True, help='Home directory where csvs of elabs live.')
    parser.add_argument('-t', required=True, help='Directory of templates to use. Name must contain fragalysis ID.')
    parser.add_argument('-i', required=True, help='SDF of all fragment hits for target.')
    parser.add_argument('--step', required=True, help='Step identifier to place elabs for. Must be present in csv name.')
    parser.add_argument('-o', required=False, help='Output csv to save logging results.')
    parser.add_argument('-p', required=False, help='Prefix for fragment names in sdf.')
    parser.add_argument('--n_cores', required=False, default=1, help='Number of cores to use.')
    parser.add_argument('-l', '--log_path', required=False, help='Path to the log directory.')
    parser.add_argument('--cutoff', action='store_true', required=False, help='Cutoff of 10,000 for placing elabs.')
    parser.add_argument('--all_csv', required=False, help='CSV of all elabs for elaboration campaign.')
    parser.add_argument('--wictor', action='store_true', required=False, help='Use if running on with Wictor in Fragmenstein')
    return parser

def shorten_elabs_csv(elabs_csv, len):
    """Shortens the elabs csv to 10,000 rows."""
    df = pd.read_csv(elabs_csv, encoding='ISO-8859-1')
    df = df[:10000]
    new_path = elabs_csv.replace(f'{len}.csv','') + "10000.csv"
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

def find_template_pdb(template_dir, frag_name):
    """Finds the template pdb with the given frag name."""
    for root, dirs, files in os.walk(template_dir):
        for file in files:
            if frag_name in file:
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

def find_elabs_csv(root, directory, frag1, frag2, step_identifier, sdf_prefix=None, add_hit_names=False):
    """Finds the elabs csv for the given directory. Adds column 'hit_names' if needed."""
    for filename in os.listdir(os.path.join(root, directory)):
        if filename.endswith('.csv') and step_identifier in filename and not filename.startswith('.'):
            if add_hit_names:
                csv_path = os.path.join(root, directory, filename)
                print(csv_path)
                df = pd.read_csv(csv_path, encoding='ISO-8859-1')
                df['hit_names'] = f"{sdf_prefix}{frag1} {sdf_prefix}{frag2}"
                df.to_csv(csv_path, index=False)
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

    for root, dirs, files in os.walk(kwargs['d']):
        for directory in dirs:
            print('DIRECTORY', directory)
            if directory == 'output':
                exit()
            done = False
            if "xtra_results" in directory or "logs" in directory:
                exit()
            if "," in directory:
                continue
            for sub_dir in os.listdir(os.path.join(root, directory)): # checks if output folder is already there, skip
                if 'output' in sub_dir:
                    done=True
            if done:
                continue
            print(f"DIRECTORY: {directory}")
            frag1 = directory.split("_")[2]+"_"+directory.split("_")[3]
            frag2 = directory.split("_")[4]+"_"+directory.split("_")[5]
            cmpd_catalog = directory.split("_")[0]+"_"+directory.split("_")[1]
            # frags_sdf is no longer needed to could be helpful for viewing even though it is output by fragmenstien,
            # I might actually remove it
            #frags_sdf = find_frags_sdf(sdf_content, root, directory, cmpd_catalog, frag1, frag2, sdf_prefix=kwargs['p'])
            frags_sdf = sdf_file_path
            print(frags_sdf)
            template_pdb = find_template_pdb(kwargs['t'], frag1)
            print(template_pdb)
            elabs_csv, len = find_elabs_csv(root, directory, frag1, frag2, step_identifier=kwargs['step'],
                                            sdf_prefix=kwargs['p'], add_hit_names=True)
            print(elabs_csv)
            if frags_sdf is None:
                print(f"Frags sdf not found for {directory}.")
            if template_pdb is None:
                print(f"Template pdb not found for {directory}.")
            if elabs_csv is None:
                print(f"Elabs csv not found for {directory}.")
            if frags_sdf and template_pdb and elabs_csv:
                os.chdir(os.path.join(root, directory)) # change to directory
                try:
                    if kwargs['cutoff'] and len > 10000:
                        print(f"CUTTING {directory} because it has more than 10,000 elabs.")
                        elabs_csv = shorten_elabs_csv(elabs_csv, len)
                    print(f"PLACING {directory}.")
                    if kwargs['wictor']:
                        subprocess.run(
                            ["fragmenstein", "laboratory", "place", "--input", frags_sdf, "--template", template_pdb,
                             "--in-table", elabs_csv, "--output",
                             os.path.join(root, directory, f"{cmpd_catalog}_{frag1}_{frag2}_output.csv"),
                             "--cores", str(n_cores), "--verbose", "--victor", "Wictor"])
                    else:
                        subprocess.run(
                            ["fragmenstein", "laboratory", "place", "--input", frags_sdf, "--template", template_pdb,
                             "--in-table", elabs_csv, "--output",
                             os.path.join(root, directory, f"{cmpd_catalog}_{frag1}_{frag2}_output.csv"),
                             "--cores", str(n_cores), "--verbose"])
                    for filename in os.listdir(os.path.join(root, directory)):
                        if filename.endswith('output.csv') and not filename.startswith('.'):
                            df = pd.read_csv(os.path.join(root, directory, filename), encoding='ISO-8859-1', index_col=0)
                            num_accept_elabs = df[df['outcome']=='acceptable'].shape[0]
                            df_output = pd.read_csv(os.path.join(root, directory, filename), encoding='ISO-8859-1', index_col=0)
                            print(f"Number of acceptable elabs: {num_accept_elabs}")
                            break
                    if kwargs['all_csv']:
                        print(f"APPENDING {directory} to all elabs csv.")
                        all_csv = kwargs['all_csv']
                        df_all = pd.read_csv(all_csv, encoding='ISO-8859-1')
                        df_all = update_placed_acceptable(df_all, directory, num_accept_elabs)
                        # Save the modified DataFrame back to the CSV
                        df_all.to_csv(all_csv, index=False)
                    print('MERGING FRAGMENSTEIN OUTPUT WITH ORIGINAL ELABS CSV.')
                    df = pd.read_csv(elabs_csv, encoding='ISO-8859-1')
                    df = df.merge(df_output, on='smiles', how='left')
                    df.to_csv(elabs_csv, index=False)
                    print("DELETING EXTRA FILES.")
                    output_dir = os.path.join(root, directory, 'output')
                    if not os.path.exists(output_dir):
                        print(f"Output directory does not exist for {directory}.")
                        continue
                    os.chdir(output_dir)  # change to output directory
                    suffixes_to_keep = ['.minimised.json', '.holo_minimised.pdb', '.minimised.mol']
                    for root1, dirs, files in os.walk(output_dir):
                        for directory2 in dirs:
                            directory2_path = os.path.join(output_dir, directory2)
                            if not os.path.exists(directory2_path):
                                print('THIS PATH DOES NOT EXIST', directory2_path)
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
                    print(f"Error placing elabs for {directory}.")
                    print(e)
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
