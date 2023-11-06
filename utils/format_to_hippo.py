import argparse
import os
import csv
import json
from typing import Dict, Any, Tuple
import glob2
import pandas as pd
import shutil

def config_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', required=True, help='Home directory for all the placements')
    parser.add_argument('-e', '--elab', required=True, help='Suffix of the elab CSV file to merge')
    parser.add_argument('-o', '--output', required=True, help='Suffix of the output CSV file to merge')
    return parser
#-------------------------#

def find_success_dirs(dir_path: str, elab: str, output: str) -> pd.DataFrame:
    """
    Find successful dirs from the merge of original elab csv and output csv from Fragmenstein.

    Returns a list of dictionaries representing the rows in the merged CSV.

    INPUT:
        dir_path: directory path of the base compound
        elab: suffix of the elab csv file
        output: suffix of the output csv file
    """
    merged_data = []
    elab_path = ''
    output_path = ''
    csv_files = glob2.glob(dir_path)
    for file in csv_files:
        if elab in file:
            elab_path = file
        elif output in file:
            output_path = file
    if elab_path == '' or output_path == '':
        return None
    elabdf = pd.read_csv(elab_path)
    outputdf = pd.read_csv(output_path)
    merged_data = elabdf.merge(outputdf, on='smiles', how='left', suffixes=('_elab', '_output'))
    #Drop rows where there is no output or an error
    success_data = merged_data.drop(merged_data[(merged_data['output'] != 'acceptable') | (merged_data['output'] != 'too moved'].index))
    return success_data

def create_success_directories(root: str, dir_path: str, merged_data: pd.DataFrame) -> Tuple[Dict[str, int], int]:
    """
    Create success directories and return a mapping of paths to number of files.
    """
    success_dict = {}
    success_dir_path = os.path.join(root, dir_path, 'success')
    if not os.path.exists(success_dir_path):
        os.makedirs(success_dir_path)
    # move dirs in success.csv to success dir
    merged_data.apply(lambda row: shutil.move(os.path.join(root, dir_path, row['name']), success_dir_path), axis=1)
    num_dirs = len(os.listdir(success_dir_path))
    success_dict[success_dir_path] = num_dirs
    merged_data.to_csv(os.path.join(success_dir_path, 'success.csv'))
    return success_dict, num_dirs

#----------------------------#
def main():
    parser = config_parser()
    args = parser.parse_args()
    total_success_dirs = []
    num_success_dirs = 0
    num_success_elabs = 0
    for root, directory, files in os.walk(args.directory): # walk through home directory, where the directories of the elab compounds are
        # Merge csvs of the original csv and output csv from Fragmenstein
        success_data = find_success_dirs(directory, args.elab, args.output)
        if success_data is None:
            print(f"{directory} DOES NOT CONTAIN THE ELAB AND OUTPUT CSV FILES")
            continue # go onto next directory of elab compounds
        # Create success dirs
        success_dict, num_dirs = create_success_directories(root, args.directory, merged_data)
        num_success_dirs += 1
        num_success_elabs += len(num_dirs)
        # Store success dirs paths and number of successes within each
        total_success_dirs.append(success_dict)
    # Save results
    with open('success_dirs.json', 'w') as f:
        json.dump(total_success_dirs, f)

if __name__ == "__main__":
    main()
