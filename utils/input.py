import ast

import pandas as pd
from ast import literal_eval

from reactionSmarts import REACTANT1_DICT, REACTANT2_DICT, reaction_smarts
from chemUtils.substructure import find_attachmentIdxs_fromMolAndSubstruct
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import FingerprintSimilarity
from constants import REACTIONS_NAMES


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

def findstep2(reaction_name, final_product, step1_product, reactants):
    """
    Given a reaction, final product, product from step before, original reactants. Return the exact reactant match to reactant 1 or 2.

    INPUT:
        reaction_name: string of reaction name
        product: string of product SMILES
        reactants: list of reactants SMILES
    OUTPUT:
        exact_reactant1: string of reactant1 SMILES matched to reactant1 smarts
        exact_reactant2: string of reactant2 SMILES matched to reactant2 smarts
    """
    if reaction_name not in REACTIONS_NAMES:
        print(f"Do not have SMARTS for this reaction: {reaction_name}\n "
              f"Please provide the SMARTS.\n"
              f"Skipping {final_product}...\n")
        return NotImplementedError
    else:
        try:
            pred_prods = {}
            smarts_pattern = reaction_smarts[reaction_name]
            reaction = AllChem.ReactionFromSmarts(smarts_pattern)
            reactant1 = step1_product # IMPORTANT
            reactant2 = [x for x in reactants if x != step1_product][0] # IMPORTANT
            if reactant2 is None:
                print('NO MATCH FOUND FOR REACTANT 2')
                return NotImplementedError
            # Find all products from both combinations
            pred_prods['r1_r2'] = [x[0] for x in reaction.RunReactants((reactant1, reactant2)) if x is not None]
            pred_prods['r2_r1'] = [x[0] for x in reaction.RunReactants((reactant2, reactant1)) if x is not None]
            # Check tanimoto
            for key, value in pred_prods.items():
                if len(value) == 0:
                    continue
                else:
                    for prod in value:
                        final_prod = Chem.MolFromSmiles(final_product)
                        Chem.SanitizeMol(final_prod)
                        Chem.SanitizeMol(prod)
                        prod_fp = AllChem.GetMorganFingerprintAsBitVect(prod, radius=2, nBits=1024)
                        final_prod_fp = AllChem.GetMorganFingerprintAsBitVect(final_prod, radius=2, nBits=1024)
                        if FingerprintSimilarity(prod_fp, final_prod_fp) == 1:
                            if key=='r1_r2':
                                exact_reactant1 = step1_product
                                exact_reactant2 = reactant2
                            if key=='r2_r1':
                                exact_reactant1 = reactant2
                                exact_reactant2 = step1_product
                            break
                print('NO MATCH FOUND FOR THIS FINAL:', final_product)
                print('REACTANTS: ', reactants)
        except Exception as e:
            print(e)
    return exact_reactant1, exact_reactant2

def editstep1(df):
    """
    For routes longer than 1 step we have to find the product of the first step and use that as the reactant for
    the next. Make whole new dataframe.

    INPUT:
        df: dataframe with columns SMILES, Reaction_name, Reactants

    OUTPUT:
        step1_output_df: dataframe with columns 'smiles', 'num_steps', 'reactants', 'rxn_order_first_to_last', 'reactants', 'dir_name'
        step2_output_df: dataframe with columns 'smiles', 'num_steps', 'reactants', 'rxn_order_first_to_last', 'reactants', 'dir_name'
    """
    step1_output_df = pd.DataFrame(columns=['smiles', 'num_steps', 'reactants', 'rxn_order_first_to_last', 'dir_name'])
    step2_output_df = pd.DataFrame(columns=['smiles', 'num_steps', 'reactants', 'rxn_order_first_to_last', 'dir_name'])
    for index, row in df.iterrows():
        reaction_name = ast.literal_eval(row['rxn_order_first_to_last'])[0].replace(' ', '_')
        reaction_name_step2 = ast.literal_eval(row['rxn_order_first_to_last'])[1].replace(' ', '_')
        if reaction_name not in REACTIONS_NAMES:
            print(f"Do not have SMARTS for this reaction: {reaction_name}\n "
                  f"Please provide the SMARTS.\n"
                  f"Skipping {row['smiles']}...\n")
            continue
        if reaction_name_step2 not in REACTIONS_NAMES:
            print(f"Do not have SMARTS for this reaction: {reaction_name_step2}\n "
                  f"Please provide the SMARTS.\n"
                  f"Skipping {row['smiles']}...\n")
            continue
        try:
            final_product = row['smiles']
            pred_prods = {}
            reactants = ast.literal_eval(row['reactants'])
            reactant1 = Chem.MolFromSmiles(reactants[0][0])
            reactant2 = Chem.MolFromSmiles(reactants[0][1])
            poss_prods = reactants[1]
            smarts_pattern = reaction_smarts[reaction_name]
            reaction = AllChem.ReactionFromSmarts(smarts_pattern)
            # Find all products from both combinations
            pred_prods['r1_r2'] = [x[0] for x in reaction.RunReactants((reactant1, reactant2)) if x is not None]
            pred_prods['r2_r1'] = [x[0] for x in reaction.RunReactants((reactant2, reactant1)) if x is not None]
            # Check tanimoto
            for key, value in pred_prods.items():
                if len(value) == 0:
                    continue
                else:
                    for prod in value:
                        for poss_prod in poss_prods:
                            poss_prod_smiles = poss_prod
                            poss_prod = Chem.MolFromSmiles(poss_prod)
                            Chem.SanitizeMol(poss_prod)
                            Chem.SanitizeMol(prod)
                            prod_fp = AllChem.GetMorganFingerprintAsBitVect(prod, radius=2, nBits=1024)
                            poss_prod_fp = AllChem.GetMorganFingerprintAsBitVect(poss_prod, radius=2, nBits=1024)
                            if FingerprintSimilarity(prod_fp, poss_prod_fp) == 1:
                                step1_output_df.loc[index, 'smiles'] = poss_prod_smiles
                                step1_output_df.loc[index, 'num_steps'] = 1
                                step1_output_df.loc[index, 'rxn_order_first_to_last'] = [ast.literal_eval(row['rxn_order_first_to_last'])[0]]
                                step1_output_df.loc[index, 'reactants'] = [reactants[0]]
                                step1_output_df.loc[index, 'dir_name'] = row['dir_name']
                                break
                print('NO MATCH FOUND FOR INDEX: ', index)
                print('POSS_PROD: ', poss_prod_smiles)
                exact_reactant1, exact_reactant2 = findstep2(reaction_name_step2, final_product, poss_prod_smiles, reactants[1])
                step2_output_df.loc[index, 'smiles'] = final_product
                step2_output_df.loc[index, 'num_steps'] = 2
                step2_output_df.loc[index, 'rxn_order_first_to_last'] = [ast.literal_eval(row['rxn_order_first_to_last'])[1]]
                step2_output_df.loc[index, 'reactant1_exact'] = exact_reactant1
                step2_output_df.loc[index, 'reactant2_exact'] = exact_reactant2
        except Exception as e:
            print(e)

    return step1_output_df

