import argparse
import os
import csv
import json
from typing import Dict, Any, Tuple
import glob2
import pandas as pd
import shutil
import unicodedata
from pathlib import Path

def config_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', required=True, help='Home directory for all the placements')
    parser.add_argument('-e', '--elab', required=True, help='Identifier of the elab CSV file to merge')
    parser.add_argument('-o', '--output', required=True, help='Identifier of the output CSV file to merge')
    parser.add_argument('--rmsd', required=False, help='RMSD threshold to accept placement. All placements '
                                                       'below this will be accepted. Default is 1.0.')
    return parser
#-------------------------#

def get_merged_data(root_path: str, dir_path: str, elab_identifier: str, output_identifier: str, rmsd_thresh: float) -> pd.DataFrame:
    """
    Gets dataframe of successful placements based on rmsd threshold.

    INPUT:
        dir_path: directory path of the base compound
        elab_suffix: suffix of the elab csv file
        output_suffix: suffix of the output csv file
    OUTPUT:
        acceptable_data: DataFrame of the acceptable placements and its metadata
        elab_path: path to the csv of elaborations
    """
    merged_data = []
    print('elab_identifier', elab_identifier)
    print('output_identifier', output_identifier)
    elab_path = ''
    output_path = ''
    full_path = os.path.join(root_path, dir_path)
    print('full_path', full_path)
    # Find all csv files in the directory
    csv_files = glob2.glob(os.path.join(full_path, '*.csv'))
    print('csv_files', csv_files)
    # Find the elab and output files based on the unique identifiers
    for file in csv_files:
        if elab_identifier in file and output_identifier not in file and 'success_moved' not in file:
            elab_path = file
        elif output_identifier in file:
            output_path = file
    print('elab_path', elab_path)
    print('output_path', output_path)
    if elab_path == '' or output_path == '':
        print(f"Could not find elab or output csv file for {dir_path}")
        return None
    elabdf = pd.read_csv(elab_path)
    elabdf.drop_duplicates(inplace=True) # always drop duplicates
    # Check if elabdf is already merged by if it contains 'outcome' column
    if 'outcome' not in elabdf.columns: # Merge
        outputdf = pd.read_csv(output_path, index_col=0)
        elabdf = elabdf.merge(outputdf, on='smiles', how='left', suffixes=('_elab', '_output'))
    elabdf.to_csv(elab_path, index=False)
    assert 'comRMSD' in elabdf.columns, f"RMSD column does not exist! Please check {elab_path}"
    # Find delta delta G column
    #ddg_identifier = 'ΔΔG'
    ddg_identifier = unicodedata.normalize('NFKC', '∆∆G')
    ddg_column = None
    for col in elabdf.columns:
        col = unicodedata.normalize('NFKC', col)
        if ddg_identifier in col:
            ddg_column = col
    if ddg_column is None:
        print(f'∆∆G was not found! Please check {elab_path}...')
    acceptable_data = elabdf[(elabdf['comRMSD'] <= rmsd_thresh) & (elabdf[ddg_column] < 0)]
    return acceptable_data, elab_path


def move_folder(row, output_path, success_dir_path, df):
    # Check if 'moved_to_success_dir' column exists in df, if not, add it
    if 'moved_to_success_dir' not in df.columns:
        df['moved_to_success_dir'] = False
    # Define the source and destination paths
    source = os.path.join(output_path, row['name'].replace('_', '-'))
    destination = os.path.join(success_dir_path, row['name'].replace('_', '-'))
    # Check if the destination directory already exists
    if not os.path.exists(destination):
        # If it does not exist, attempt to move the folder
        try:
            shutil.move(source, destination)
            # Update the DataFrame if move is successful
            df.at[row.name, 'moved_to_success_dir'] = True
        except Exception as e:
            # If there's an error, print it, but continue with the next iteration
            print(f"Could not move {source} to {destination}: {e}")
    else:
        # If it exists, print a message and continue
        df.at[row.name, 'moved_to_success_dir'] = True
        print(f"Directory {destination} already exists, skipping.")

def create_success_directories(root: str, dir: str, acceptable_data: pd.DataFrame, elab_path: str) -> Tuple[Dict[str, int], int]:
    """
    Create success directories and return a mapping of paths to number of files. Will do an exhaustive check if success
    dirs already exist and will move the ones that have not been created yet.

    INPUT:
        root: root path to parent directory
        dir: directory name
        acceptable_data: dataframe of the metadata of the acceptable placements
        elab_path: path to the elaboration csv
    """
    success_dict = {}
    success_dir_path = os.path.join(root, dir, 'success')
    if not os.path.exists(success_dir_path):
        os.makedirs(success_dir_path)
    output_path = os.path.join(root, dir) + '/output'
    # move dirs in success.csv to success dir
    try:
        acceptable_data.apply(lambda row: move_folder(row, output_path=output_path, success_dir_path=success_dir_path, df=acceptable_data), axis=1)
    except Exception as e:
        print(e)
    num_dirs = len(os.listdir(success_dir_path))
    success_dict[success_dir_path] = num_dirs
    acceptable_data.to_csv(os.path.join(success_dir_path, 'success.csv'))
    # Merge elabdf with acceptable data
    elabdf = pd.read_csv(elab_path)
    acceptable_data_subset = acceptable_data[['smiles', 'moved_to_success_dir']]
    elab_path2 = Path(elab_path)
    success_path_identifier = '_success_moved'
    # Check if 'moved_to_success_dir' already exists, if so, delete. Also delete suffix of
    # success_moved since adding it again will be redundant.
    if 'moved_to_success_dir' in elabdf.columns:
        del elabdf['moved_to_success_dir']
        # Remove the identifier from the stem
        new_stem = elab_path2.stem.replace(success_path_identifier, '')
        elab_path2 = elab_path2.with_name(new_stem + elab_path2.suffix)
    df = elabdf.merge(acceptable_data_subset, on='smiles', how='left', suffixes=('_elab', '_output'))
    df.drop_duplicates(inplace=True)
    success_dir_merge_path = elab_path2.with_stem(f"{elab_path2.stem}{success_path_identifier}")
    # Save merged data to new csv
    df.to_csv(success_dir_merge_path, index=False)
    # Save elab dataframe to same path
    elabdf.to_csv(elab_path, index=False)
    return success_dict, num_dirs

#----------------------------#
def main():
    parser = config_parser()
    args = parser.parse_args()
    total_success_dirs = []
    num_success_dirs = 0
    num_success_elabs = 0
    for root, directories, files in os.walk(args.directory): # walk through home directory, where the directories of the elab compounds are
        # Only interested in the top-level directories within args.directory
        if root != args.directory:
            # If the root is not the directory we started with, it's a subdirectory, so skip it
            continue
        for directory in directories:
            if directory == 'success' or directory == 'output' or 'xtra_results' in directory:
                continue
            # Merge csvs of the original csv and output csv from Fragmenstein (if not already done)
            try:
                if args.rmsd is None:
                    rmsd_thresh = 1.0
                else:
                    rmsd_thresh = float(args.rmsd)
                acceptable_data, elab_path = get_merged_data(root, directory, args.elab, args.output, rmsd_thresh)
                if acceptable_data is None:
                    print(f"{directory} DOES NOT CONTAIN THE ELAB AND OUTPUT CSV FILES")
                    continue # go onto next directory of elab compounds
                # Create success dirs
                success_dict, num_dirs = create_success_directories(root, directory, acceptable_data, elab_path)
                num_success_dirs += 1
                num_success_elabs += num_dirs
                # Store success dirs paths and number of successes within each
                total_success_dirs.append(success_dict)
                print(f"SUCCESS: {directory}")
            except Exception as e:
                print(f"Error in {directory}: {e}")
        break
    # Save results
    with open((os.path.join(args.directory, 'success_dirs.json')), 'w') as f:
        json.dump(total_success_dirs, f)

if __name__ == "__main__":
    main()
