#!/usr/bin/env python3
"""
syndirella.justretroquery.py

This module provides functions to output retrosynthesis queries for a given list of scaffolds.
"""
import json
import logging
import os.path
from typing import Dict, List, Tuple, Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

from syndirella.database.Postera import Postera
from syndirella.utils.error import APIQueryError, NoSynthesisRoute
from syndirella.constants import RetrosynthesisTool, DEFAULT_RETROSYNTHESIS_TOOL
from .cli_defaults import cli_default_settings

logger = logging.getLogger(__name__)
with open(cli_default_settings['rxn_smarts_path']) as f:
    reaction_smarts = json.load(f)
reaction_smarts_names: List[str] = list(reaction_smarts.keys())

# Load SMIRKS patterns for validation
reaction_smirks_dict = {name: data['smirks'] for name, data in reaction_smarts.items()}


def save_df(df: pd.DataFrame, output_dir: str, csv_path: str, retro_tool: str = None) -> Tuple[str, str]:
    """
    Save the DataFrame to the output directory as CSV and as gzip-compressed pickle.
    
    Args:
        df: DataFrame to save
        output_dir: Output directory path
        csv_path: Original CSV path (used for naming)
        retro_tool: Optional retrosynthesis tool name to include in filename
    
    Returns:
        Tuple of (path to saved CSV file, path to saved pickle file)
    """
    csv_basename = os.path.basename(csv_path)
    base_name = csv_basename.replace('.csv', '')
    
    # Create filename prefix with tool name if provided
    if retro_tool:
        prefix = f'justretroquery_{retro_tool.lower()}_'
    else:
        prefix = 'justretroquery_'
    
    # Save CSV format
    csv_output_name = f'{prefix}{base_name}.csv'
    csv_output_path = os.path.join(output_dir, csv_output_name)
    df.to_csv(csv_output_path, index=False)
    
    # Save gzip-compressed pickle format
    pkl_output_name = f'{prefix}{base_name}.pkl.gz'
    pkl_output_path = os.path.join(output_dir, pkl_output_name)
    df.to_pickle(pkl_output_path)
    
    logger.info(f"Saved DataFrame to CSV: {csv_output_path}")
    logger.info(f"Saved DataFrame to pickle: {pkl_output_path}")
    
    return csv_output_path, pkl_output_path


def validate_reaction_with_smirks(reaction_name: str, reactant_smiles: Tuple[str, ...], 
                                   expected_product_smiles: str) -> Tuple[bool, Optional[str]]:
    """
    Validate that applying a Syndirella SMIRKS reaction to reactants produces the expected product.
    
    Args:
        reaction_name: Name of the reaction (must be in Syndirella's reaction library)
        reactant_smiles: Tuple of reactant SMILES strings
        expected_product_smiles: Expected product SMILES from the retrosynthesis route
    
    Returns:
        Tuple of (is_valid, actual_product_smiles or None)
        - is_valid: True if the reaction produces a product matching expected_product_smiles (by InChI-key)
        - actual_product_smiles: SMILES of the product produced, or None if no product or validation failed
    """
    # Check if reaction is in Syndirella's library
    if reaction_name not in reaction_smirks_dict:
        logger.debug(f"Reaction {reaction_name} not in Syndirella library, skipping validation")
        return False, None
    
    try:
        # Get SMIRKS pattern
        smirks = reaction_smirks_dict[reaction_name]
        
        # Create reaction from SMIRKS
        reaction = ReactionFromSmarts(smirks)
        
        # Convert reactant SMILES to molecules
        reactant_mols = []
        for smiles in reactant_smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.debug(f"Failed to parse reactant SMILES: {smiles}")
                return False, None
            reactant_mols.append(mol)
        
        # Apply reaction
        products = reaction.RunReactants(tuple(reactant_mols))
        
        if not products or len(products) == 0:
            logger.debug(f"No products produced for reaction {reaction_name}")
            return False, None
        
        # Get expected product InChI-key
        expected_mol = Chem.MolFromSmiles(expected_product_smiles)
        if expected_mol is None:
            logger.debug(f"Failed to parse expected product SMILES: {expected_product_smiles}")
            return False, None
        expected_inchi_key = Chem.MolToInchiKey(expected_mol)
        
        # Check if any product matches the expected product
        for product_tuple in products:
            for product_mol in product_tuple:
                try:
                    # Sanitize product
                    Chem.SanitizeMol(product_mol)
                    product_inchi_key = Chem.MolToInchiKey(product_mol)
                    
                    if product_inchi_key == expected_inchi_key:
                        product_smiles = Chem.MolToSmiles(product_mol)
                        return True, product_smiles
                except Exception as e:
                    logger.debug(f"Error sanitizing/processing product: {e}")
                    continue
        
        # If we get here, no product matched
        if products:
            # Return the first product SMILES even if it doesn't match (for debugging)
            try:
                first_product = products[0][0]
                Chem.SanitizeMol(first_product)
                actual_product_smiles = Chem.MolToSmiles(first_product)
                logger.debug(f"Product mismatch for {reaction_name}: expected {expected_product_smiles}, got {actual_product_smiles}")
                return False, actual_product_smiles
            except:
                pass
        
        return False, None
        
    except Exception as e:
        logger.debug(f"Error validating reaction {reaction_name}: {e}")
        return False, None


def format_routes(routes: List[Dict[str, List[Dict[str, str]]]]) -> Dict:
    """
    Gets the top 5 passing routes from the routes. Formats them into a dictionary with routes names as keys, also adds
    other field of routeX_names.
    
    Args:
        routes: List of route dictionaries, each with a 'reactions' key containing list of reaction dicts.
                Routes should be pre-sorted by score (best first).
    
    Returns:
        Dictionary with route information (route0, route0_names, route0_CAR, etc. for up to 5 routes)
    """
    passing_routes = {}
    
    # Filter to only routes with reactions
    valid_routes = [route for route in routes if route.get('reactions') and len(route['reactions']) > 0]
    
    if len(valid_routes) == 0:
        logger.error("No routes with reactions found.")
        return passing_routes
    
    # Take top 5 routes (they should already be sorted by score)
    top_routes = valid_routes[:5]
    
    for i, route in enumerate(top_routes):
        # get the reactions
        reactions: List[Dict[str, str]] = route['reactions']
        # in original reactions list, replace spaces with underscores
        for reaction in reactions:
            reaction['name'] = reaction['name'].replace(" ", "_")
        
        # Validate reactions with Syndirella SMIRKS
        validated_reactions = []
        validation_results = []
        for reaction in reactions:
            reaction_name = reaction['name']
            reactant_smiles = reaction.get('reactantSmiles', ())
            expected_product_smiles = reaction.get('productSmiles', '')
            
            # Validate if reaction is in Syndirella library and has required data
            is_valid = None
            actual_product = None
            if reaction_name in reaction_smarts_names and reactant_smiles and expected_product_smiles:
                # Convert reactantSmiles to tuple if it's not already
                if isinstance(reactant_smiles, tuple):
                    reactant_tuple = reactant_smiles
                elif isinstance(reactant_smiles, list):
                    reactant_tuple = tuple(reactant_smiles)
                else:
                    reactant_tuple = (reactant_smiles,)
                
                is_valid, actual_product = validate_reaction_with_smirks(
                    reaction_name, reactant_tuple, expected_product_smiles
                )
            
            # Add validation info to reaction dict
            reaction['smirks_validated'] = is_valid
            if actual_product:
                reaction['actual_product_smiles'] = actual_product
            validated_reactions.append(reaction)
            validation_results.append(is_valid)
        
        # Use route index (0-4) for naming
        passing_routes[f'route{i}'] = validated_reactions
        passing_routes[f'route{i}_names']: List[str] = [reaction['name'] for reaction in validated_reactions]
        passing_routes[f'route{i}_CAR']: bool = all(
            [True if reaction['name'] in reaction_smarts_names else False for reaction in validated_reactions])
        passing_routes[f'route{i}_non_CAR']: List[str] | None = ([reaction['name'] for reaction in validated_reactions if
                                                                   reaction['name'] not in reaction_smarts_names]
                                                                  or None)
        # Add validation summary
        passing_routes[f'route{i}_smirks_validated'] = all([v for v in validation_results if v is not None]) if any(v is not None for v in validation_results) else None
        passing_routes[f'route{i}_num_validated'] = sum([1 for v in validation_results if v is True])
        passing_routes[f'route{i}_num_failed_validation'] = sum([1 for v in validation_results if v is False])
        
        # Add CAR_and_validated: True only if route is CAR AND all reactions passed validation
        car_status = passing_routes[f'route{i}_CAR']
        validation_status = passing_routes[f'route{i}_smirks_validated']
        passing_routes[f'route{i}_CAR_and_validated'] = (car_status is True and validation_status is True)
        
        logger.info(f"Route {i} has {len(validated_reactions)} reaction(s)")
        logger.info(f"Route {i} CAR: {passing_routes[f'route{i}_CAR']}")
        logger.info(f"Route {i} non-CAR: {passing_routes[f'route{i}_non_CAR']}")
        if passing_routes[f'route{i}_smirks_validated'] is not None:
            logger.info(f"Route {i} SMIRKS validation: {passing_routes[f'route{i}_num_validated']}/{len([v for v in validation_results if v is not None])} reactions validated successfully")
        logger.info(f"Route {i} CAR_and_validated: {passing_routes[f'route{i}_CAR_and_validated']}")
    
    # Calculate number of routes (count keys that are just 'routeX' without suffixes)
    num_routes = len([k for k in passing_routes.keys() if k.startswith('route') and k.replace('route', '').isdigit()])
    logger.info(f"Formatted {num_routes} routes (top {len(top_routes)} scored routes)")
    return passing_routes


def retro_search(scaffold: str, retro_tool: RetrosynthesisTool = DEFAULT_RETROSYNTHESIS_TOOL,
                 aizynth_search=None, postera=None) -> pd.DataFrame:
    """
    Perform retrosynthesis search on the given scaffold and formats outputs.
    
    Args:
        scaffold: SMILES string of the scaffold to search
        retro_tool: Retrosynthesis tool to use (Manifold or AiZynthFinder)
        aizynth_search: Optional pre-initialized AiZynthManager instance (for efficiency)
        postera: Optional pre-initialized Postera instance (for efficiency)
    
    Returns:
        DataFrame with retrosynthesis route information
    """
    logger.info(f"Performing retrosynthesis search for {scaffold} using {retro_tool}")
    
    if retro_tool == RetrosynthesisTool.MANIFOLD:
        # Use Postera/Manifold for retrosynthesis
        if postera is None:
            postera = Postera()
        routes: List[Dict[str, List[Dict[str, str]]]] | None = postera.perform_route_search(scaffold)
        if routes is None:
            logger.critical(f"API retrosynthesis query failed for {scaffold}.")
            raise APIQueryError(message=f"API retrosynthesis query failed for {scaffold}.", smiles=scaffold)
    elif retro_tool == RetrosynthesisTool.AIZYNTHFINDER:
        # Use AiZynthFinder for retrosynthesis
        if aizynth_search is None:
            try:
                from syndirella.aizynth.AiZynthManager import AiZynthManager
                aizynth_search = AiZynthManager()
            except ImportError:
                logger.error("AiZynthFinder not available. Falling back to Manifold.")
                if postera is None:
                    postera = Postera()
                routes: List[Dict[str, List[Dict[str, str]]]] | None = postera.perform_route_search(scaffold)
                if routes is None:
                    logger.critical(f"API retrosynthesis query failed for {scaffold}.")
                    raise APIQueryError(message=f"API retrosynthesis query failed for {scaffold}.", smiles=scaffold)
                formatted_routes: Dict = format_routes(routes)
                logger.info(f"Found {len(formatted_routes)} formatted routes")
                to_add = {'smiles': scaffold}
                for route, details in formatted_routes.items():
                    to_add[route] = details
                return pd.DataFrame([to_add])
        
        routes: List[Dict[str, List[Dict[str, str]]]] = aizynth_search.perform_route_search(
            target_smiles=scaffold,
            matching_strategy='best_overall',
            max_routes_per_cluster=5  # Get up to 5 routes per cluster to ensure we have enough for top 5 overall
        )
        # Convert AiZynthFinder format to expected format
        routes = _convert_aizynth_routes(routes)
        # Limit to top 5 routes (they're already sorted by score from get_best_routes)
        routes = routes[:5]
    else:
        raise ValueError(f"Unsupported retrosynthesis tool: {retro_tool}")
    
    formatted_routes: Dict = format_routes(routes)
    logger.info(f"Found {len(formatted_routes)} formatted routes")
    to_add = {'smiles': scaffold}
    for route, details in formatted_routes.items():
        to_add[route] = details
    return pd.DataFrame([to_add])  # each dictionary is a single row


def _convert_aizynth_routes(aizynth_routes: List[Dict]) -> List[Dict]:
    """
    Convert AiZynthFinder route format to the expected format.
    
    Args:
        aizynth_routes: Routes from AiZynthFinder
        
    Returns:
        Converted routes in expected format
    """
    # TODO: Fix this
    converted_routes = []
    for route in aizynth_routes:
        # Convert AiZynthFinder format to expected format
        # This is a placeholder - actual conversion depends on AiZynthFinder output format
        converted_route = {
            'reactions': route.get('reactions', [])
        }
        converted_routes.append(converted_route)
    return converted_routes


def process_df(df: pd.DataFrame, retro_tool: RetrosynthesisTool = DEFAULT_RETROSYNTHESIS_TOOL,
               output_dir: str = None, csv_path: str = None):
    """
    Process the input DataFrame and create output df with retrosynthesis information.
    Saves incrementally after each scaffold is processed.
    
    Args:
        df: Input DataFrame with 'smiles' column
        retro_tool: Retrosynthesis tool to use
        output_dir: Optional output directory for incremental saves
        csv_path: Optional CSV path for incremental saves
    """
    logger.info(f"Processing DataFrame with len {len(df)} using {retro_tool}")
    
    # Initialize the retrosynthesis tool once for efficiency
    aizynth_search = None
    postera = None
    
    if retro_tool == RetrosynthesisTool.AIZYNTHFINDER:
        try:
            from syndirella.aizynth.AiZynthManager import AiZynthManager
            logger.info("Initializing AiZynthManager (this may take a moment)...")
            aizynth_search = AiZynthManager()
            logger.info("AiZynthManager initialized successfully")
        except ImportError:
            logger.warning("AiZynthFinder not available. Will fall back to Manifold for each scaffold.")
    elif retro_tool == RetrosynthesisTool.MANIFOLD:
        postera = Postera()
    
    route_infos = []
    for i, scaffold in enumerate(df['smiles']):
        logger.info(f"Processing scaffold {i+1}/{len(df)}: {scaffold}")
        try:
            route_info: pd.DataFrame = retro_search(scaffold, retro_tool, aizynth_search=aizynth_search, postera=postera)
            route_infos.append(route_info)
        except NoSynthesisRoute as e:
            error_msg = f"No synthesis route found by {retro_tool.value}"
            logger.warning(f"{error_msg} for scaffold {i+1}/{len(df)}: {scaffold}. Continuing with next scaffold.")
            # Create a minimal DataFrame with error information
            route_info = pd.DataFrame([{
                'smiles': scaffold,
                'error': error_msg,
                'error_type': 'NoSynthesisRoute'
            }])
            route_infos.append(route_info)
        except APIQueryError as e:
            error_msg = f"API retrosynthesis query failed: {str(e)}"
            logger.error(f"{error_msg} for scaffold {i+1}/{len(df)}: {scaffold}. Continuing with next scaffold.")
            # Create a minimal DataFrame with error information
            route_info = pd.DataFrame([{
                'smiles': scaffold,
                'error': error_msg,
                'error_type': 'APIQueryError'
            }])
            route_infos.append(route_info)
        except Exception as e:
            error_msg = f"Unexpected error: {str(e)}"
            logger.error(f"{error_msg} for scaffold {i+1}/{len(df)}: {scaffold}. Continuing with next scaffold.")
            # Create a minimal DataFrame with error information
            route_info = pd.DataFrame([{
                'smiles': scaffold,
                'error': error_msg,
                'error_type': type(e).__name__
            }])
            route_infos.append(route_info)
        
        # Save incrementally after each scaffold if output paths provided
        if output_dir is not None and csv_path is not None:
            try:
                route_info_df = pd.concat(route_infos)
                # Merge with full df - unprocessed scaffolds will have NaN route columns
                merged_df = df.merge(route_info_df, on='smiles', how='left')
                merged_df.reset_index(drop=True, inplace=True)
                csv_saved, pkl_saved = save_df(merged_df, output_dir, csv_path, retro_tool=retro_tool.value)
                logger.info(f"Incrementally saved results after scaffold {i+1}/{len(df)} (CSV: {csv_saved}, pickle: {pkl_saved})")
            except Exception as e:
                logger.warning(f"Failed to save incrementally after scaffold {i+1}: {e}")
    
    # format the final DataFrame
    route_info_df = pd.concat(route_infos)
    merged_df = df.merge(route_info_df, on='smiles', how='outer')
    merged_df.reset_index(drop=True, inplace=True)
    return merged_df


#######################################

def run_justretroquery(settings: Dict):
    """
    Run the justretroquery pipeline with the given settings.
    """
    # all you need is input csv with just a column of smiles
    try:
        csv_path: str = settings['input']
        output_dir: str = settings['output']
        retro_tool = RetrosynthesisTool.from_string(settings.get('retro_tool', DEFAULT_RETROSYNTHESIS_TOOL.value))
    except KeyError as e:
        raise KeyError(f"Missing critical argument to run justretroquery: {e}")

    logger.info(f"Starting justretroquery with retrosynthesis tool: {retro_tool}")
    
    df: pd.DataFrame = pd.read_csv(csv_path)

    # Process with incremental saving
    df: pd.DataFrame = process_df(df, retro_tool, output_dir=output_dir, csv_path=csv_path)

    # Final save (in case incremental saves weren't used or for final update)
    csv_saved, pkl_saved = save_df(df, output_dir, csv_path, retro_tool=retro_tool.value)

    logger.info(f"Final DataFrame saved to CSV: {csv_saved}")
    logger.info(f"Final DataFrame saved to pickle: {pkl_saved}")

    logger.info('Justretroquery execution completed successfully.')
