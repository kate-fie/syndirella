import ast

import pandas as pd
from ast import literal_eval


def validate_routes(smiles, reaction_name, reactants):
    """Function to validate the data types"""
    if not isinstance(smiles, str):
        raise ValueError(f"Expected SMILES type as string but got {type(smiles)}")

    if not isinstance(reaction_name, str):
        raise ValueError(
            f"Expected reaction_name as string but got {type(reaction_name)} or incorrect length/value type")

    if not isinstance(reactants, tuple):
        raise ValueError(
            f"Expected reactants as a tuple but got {type(reactants)}")

    return True


def process_routes(df, smiles_col_idx, reaction_name_col_idx, reactants_col_idx, results_dir, output_filename):
    """
    Processes the dataframe or csv and returns a new CSV in the specified format.
    """
    df = df.copy()
    # Check if indices are valid
    if any(idx >= df.shape[1] for idx in [smiles_col_idx, reaction_name_col_idx, reactants_col_idx]):
        raise ValueError("Invalid column index provided")

    # Extract and evaluate columns
    df['SMILES'] = df.iloc[:, smiles_col_idx]
    df['reaction_name'] = df.iloc[:, reaction_name_col_idx].apply(literal_eval).apply(lambda x: x[0])
    df['reactants'] = df.iloc[:, reactants_col_idx].apply(literal_eval).apply(lambda x: x[0])

    # Validate the data
    for _, row in df.iterrows():
        validate_routes(row['SMILES'], row['reaction_name'], row['reactants'])

    # Reorder the dataframe
    df = df[['SMILES', 'reaction_name', 'reactants']]

    # Write to CSV
    output_path = f"{results_dir}/{output_filename}"
    df.to_csv(output_path, index=False)
    print(f"Processed data saved at: {output_path}")
