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

from typing import List, Dict, Any, Optional

def config_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', required=True, help='Home directory where csvs of elabs live.')
    parser.add_argument('-t', required=True, help='Directory of templates to use. Name must contain fragalysis ID.')
    parser.add_argument('-i', required=True, help='SDF of all fragment hits for target.')
    parser.add_argument('-o', required=False, help='Output csv to save logging results.')
    parser.add_argument('-p', required=False, help='Prefix for fragment names in sdf.')
    parser.add_argument('--n_cores', required=False, default=1, help='Number of cores to use.')
    parser.add_argument('-l', '--log_path', required=False, help='Path to the log directory.')
    return parser

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

def find_elabs_csv(root, directory, frag1, frag2, sdf_prefix=None, add_hit_names=False):
    """Finds the elabs csv for the given directory. Adds column 'hit_names' if needed."""
    for filename in os.listdir(os.path.join(root, directory)):
        if filename.endswith('.csv'):
            if add_hit_names:
                csv_path = os.path.join(root, directory, filename)
                df = pd.read_csv(csv_path)
                df['hit_names'] = f"{sdf_prefix}{frag1} {sdf_prefix}{frag2}"
                df.to_csv(csv_path, index=False)
            return csv_path, len(df)
    return None, None

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
            if "xtra_results" in directory:
                exit()
            if "," in directory:
                continue
            for sub_dir in os.listdir(os.path.join(root, directory)): # checks if output folder is already there, skip
                if 'output' in sub_dir:
                    continue

            print(f"DIRECTORY: {directory}")
            frag1 = directory.split("_")[2]+"_"+directory.split("_")[3]
            frag2 = directory.split("_")[4]+"_"+directory.split("_")[5]
            cmpd_catalog = directory.split("_")[0]+"_"+directory.split("_")[1]
            # frags_sdf is no longer needed to could be helpful for viewing even though it is output by fragmenstien,
            # I might actually remove it
            #frags_sdf = find_frags_sdf(sdf_content, root, directory, cmpd_catalog, frag1, frag2, sdf_prefix=kwargs['p'])
            frags_sdf = sdf_file_path
            template_pdb = find_template_pdb(kwargs['t'], frag1)
            elabs_csv, len = find_elabs_csv(root, directory, frag1, frag2, sdf_prefix=kwargs['p'], add_hit_names=True)
            if frags_sdf is None:
                print(f"Frags sdf not found for {directory}.")
            if template_pdb is None:
                print(f"Template pdb not found for {directory}.")
            if elabs_csv is None:
                print(f"Elabs csv not found for {directory}.")
            if frags_sdf and template_pdb and elabs_csv:
                os.chdir(os.path.join(root, directory)) # change to directory
                try:
                    print(f"PLACING {directory}.")
                    subprocess.run(
                        ["fragmenstein", "laboratory", "place", "--input", frags_sdf, "--template", template_pdb,
                         "--in-table", elabs_csv, "--output",
                         os.path.join(root, directory, f"{cmpd_catalog}_{frag1}_{frag2}_output.csv"),
                         "--cores", n_cores, "--verbose"])
                except Exception as e:
                    print(f"Error placing elabs for {directory}.")
                    continue
            else:
                continue

def main():
    parser = config_parser()
    # load
    settings: Dict[str, Any] = vars(parser.parse_args())
    # Redirect stdout and stderr to log file
    # Generate a timestamp string and append to the log file path
    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    log_path_with_timestamp = f"{settings['log']}frag_place_{timestamp}.log"
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
