#!/usr/bin/env python3
"""
syndirella.run_pipeline.py

This script contains the main pipeline for syndirella.
"""
import os
from typing import List, Tuple, Dict, Any
import pandas as pd
from rdkit import Chem
import datetime
import ast

from syndirella.Cobbler import Cobbler
from syndirella.cobblers_workshop.CobblersWorkshop import CobblersWorkshop
from syndirella.slipper.Slipper import Slipper


def _format_additional_info(df: pd.DataFrame, additional_columns: List[str]) -> Dict[str, Any]:
    """
    This function is used to format the additional info from the dataframe into a dictionary.
    """
    additional_info = {}
    for col in additional_columns:
        additional_info[col] = df[col]
    return additional_info


def _assert_csv(csv_path: str) -> str:
    """
    Make sure that the csv path exists.
    """
    assert os.path.exists(csv_path), "The csv path does not exist."
    # make sure that the csv contains the correct headers
    df = pd.read_csv(csv_path)
    required_columns = ['smiles', 'compound_set', 'hits']
    for col in required_columns:
        assert col in df.columns, f"The csv must contain the column {col}."
    # TODO: could assert that the hits column is a space seperated list of strings
    return csv_path


def _assert_manual_df(df: pd.DataFrame) -> None:
    """
    Assert that the manual dataframe is in the correct format.
    """
    required_columns = ['smiles', 'reactants', 'reaction_names', 'num_steps']
    for col in required_columns:
        assert col in df.columns, f"If doing a manual route, the csv must contain the column {col}."
    # try to read the reactants and reaction_names columns
    for index, row in df.iterrows():
        reactants = ast.literal_eval(row['reactants'])
        reaction_names = ast.literal_eval(row['reaction_names'])
        num_steps = row['num_steps']
        assert len(reactants) == len(reaction_names), \
            "The reactants and reaction_names columns must have the same number of elements."
        assert type(reactants) == list, "The reactants column must be a list of strings."
        assert type(reaction_names) == list, "The reaction_names column must be a list of strings."
        assert type(num_steps) == int, "The num_steps column must be a integer."


def _elaborate_compound_with_manual_routes(product: str,
                                           reactants: List[Tuple[str]],
                                           reaction_names: List[str],
                                           num_steps: int,
                                           hits: str,
                                           template_path: str,
                                           hits_path: str,
                                           batch_num: int,
                                           output_dir: str,
                                           additional_info: Dict[str, Any] = []):
    """
    This function is used to elaborate a single compound using a manually defined route.
    """
    mol = Chem.MolFromSmiles(product)
    assert mol, f"Could not create a molecule from the smiles {product}."
    # convert hits to a list
    hits = hits.split()
    # Just creating cobblers workshop directly
    workshop = CobblersWorkshop(product=product,
                                reactants=reactants,
                                reaction_names=reaction_names,
                                num_steps=num_steps,
                                output_dir=output_dir,
                                filter=False)
    final_library = workshop.get_final_library()
    slipper = Slipper(final_library,
                      template_path,
                      hits_path,
                      hits,
                      batch_num,
                      additional_info=additional_info)
    slipper.get_products()
    slipper.place_products()
    slipper.clean_up_placements()
    print(f"Finished elaborating compound {product} at {datetime.datetime.now()}.")
    print()


def _elaborate_compound_full_auto(product: str,
                                  hits: str,
                                  template_path: str,
                                  hits_path: str,
                                  batch_num: int,
                                  output_dir: str,
                                  additional_info: Dict[str, Any] = []):
    """
    This function is used to elaborate a single compound.
    """
    mol = Chem.MolFromSmiles(product)
    assert mol, f"Could not create a molecule from the smiles {product}."
    # convert hits to a list
    hits = hits.split()
    # create the cobbler
    cobbler = Cobbler(base_compound=product,
                      output_dir=output_dir,
                      additional_info=additional_info)
    # get the cobbler workshops
    cobbler_workshops: List[CobblersWorkshop] = cobbler.get_routes()
    for workshop in cobbler_workshops:
        final_library = workshop.get_final_library()
        slipper = Slipper(final_library,
                          template_path,
                          hits_path,
                          hits,
                          batch_num)
        slipper.get_products()
        slipper.place_products()
    print(f"Finished elaborating compound {product} at {datetime.datetime.now()}.")
    print()


def run_pipeline(csv_path: str,
                 output_dir: str,
                 template_path: str,
                 hits_path: str,
                 batch_num: int,
                 additional_columns: List[str] = [],
                 manual_routes: bool = False,
                 ):
    """
    Run the whole syndirella pipeline! 👑
    """
    csv_path = _assert_csv(csv_path)
    df = pd.read_csv(csv_path)
    # assert that the df contains additional_columns
    for col in additional_columns:
        assert col in df.columns, f"The csv must contain the column {col}."
    additional_info = _format_additional_info(df, additional_columns)
    if not manual_routes:
        print("Running the full auto pipeline.")
        for index, row in df.iterrows():  # could make this a parallel for loop
            # TODO: finish this function
            _elaborate_compound_full_auto(product=row['smiles'],
                                          hits=row['hits'],
                                          template_path=template_path,
                                          hits_path=hits_path,
                                          batch_num=batch_num,
                                          output_dir=output_dir,
                                          additional_info=additional_info)
    else:
        print("Running the pipeline with manual routes.")
        _assert_manual_df(df)
        for index, row in df.iterrows():
            _elaborate_compound_with_manual_routes(product=row['smiles'],
                                                   reactants=ast.literal_eval(row['reactants']),
                                                   reaction_names=ast.literal_eval(row['reaction_names']),
                                                   num_steps=row['num_steps'],
                                                   hits=row['hits'],
                                                   template_path=template_path,
                                                   hits_path=hits_path,
                                                   batch_num=batch_num,
                                                   output_dir=output_dir,
                                                   additional_info=additional_info)

    print("Pipeline complete.")