#!/usr/bin/env python3
"""
syndirella.justretroquery.py

This module provides functions to output retrosynthesis queries for a given list of scaffolds.
"""
import json
import logging
import os.path
from typing import Dict, List

import pandas as pd

from syndirella.database.Postera import Postera
from syndirella.utils.error import APIQueryError
from syndirella.constants import RetrosynthesisTool, DEFAULT_RETROSYNTHESIS_TOOL
from .cli_defaults import cli_default_settings

logger = logging.getLogger(__name__)
with open(cli_default_settings['rxn_smarts_path']) as f:
    reaction_smarts = json.load(f)
reaction_smarts_names: List[str] = list(reaction_smarts.keys())


def save_df(df: pd.DataFrame, output_dir: str, csv_path: str) -> str:
    """
    Save the DataFrame to the output directory.
    """
    csv_basename = os.path.basename(csv_path)
    pkl_basename = csv_basename.replace('.csv', '.pkl.gz')
    saved_path = os.path.join(output_dir, f'justretroquery_{pkl_basename}')
    df.to_pickle(saved_path)
    return saved_path


def format_routes(routes: List[Dict[str, List[Dict[str, str]]]]) -> Dict:
    """
    Gets the top 5 passing routes from the routes. Formats them into a dictionary with routes names as keys, also adds
    other field of routeX_names.
    """
    passing_routes = {}
    n_of_rxns: List[int] = [i for i, route in enumerate(routes) if len(route['reactions']) > 0]
    if len(n_of_rxns) == 0:
        logger.error("No routes with reactions found.")
        return passing_routes
    min_i = min(n_of_rxns)
    for i, route in enumerate(routes):
        # get the reactions
        reactions: List[Dict[str, str]] = route['reactions']
        # check if all reactions are in the reactions
        if len(reactions) == 0:
            continue
        # adjust i to be in range(5)
        adjusted_i = i - min_i
        # in original reactions list, replace spaces with underscores
        for reaction in reactions:
            reaction['name'] = reaction['name'].replace(" ", "_")
        passing_routes[f'route{adjusted_i}'] = reactions
        passing_routes[f'route{adjusted_i}_names']: List[str] = [reaction['name'] for reaction in reactions]
        passing_routes[f'route{adjusted_i}_CAR']: bool = all(
            [True if reaction['name'] in reaction_smarts_names else False for reaction in reactions])
        passing_routes[f'route{adjusted_i}_non_CAR']: List[str] | None = ([reaction['name'] for reaction in reactions if
                                                                           reaction[
                                                                               'name'] not in reaction_smarts_names]
                                                                          or None)
        logger.info(f"Route {adjusted_i} has {len(reactions)} reaction(s)")
        logger.info(f"Route {adjusted_i} CAR: {passing_routes[f'route{adjusted_i}_CAR']}")
        logger.info(f"Route {adjusted_i} non-CAR: {passing_routes[f'route{adjusted_i}_non_CAR']}")
        if len(passing_routes) == 20:  # 5 routes with 3 fields each
            break
    return passing_routes


def retro_search(scaffold: str, retro_tool: RetrosynthesisTool = DEFAULT_RETROSYNTHESIS_TOOL) -> pd.DataFrame:
    """
    Perform retrosynthesis search on the given scaffold and formats outputs.
    
    Args:
        scaffold: SMILES string of the scaffold to search
        retro_tool: Retrosynthesis tool to use (Manifold or AiZynthFinder)
    
    Returns:
        DataFrame with retrosynthesis route information
    """
    logger.info(f"Performing retrosynthesis search for {scaffold} using {retro_tool}")
    
    if retro_tool == RetrosynthesisTool.MANIFOLD:
        # Use Postera/Manifold for retrosynthesis
        postera = Postera()
        routes: List[Dict[str, List[Dict[str, str]]]] | None = postera.perform_route_search(scaffold)
        if routes is None:
            logger.critical(f"API retrosynthesis query failed for {scaffold}.")
            raise APIQueryError(message=f"API retrosynthesis query failed for {scaffold}.", smiles=scaffold)
    elif retro_tool == RetrosynthesisTool.AIZYNTHFINDER:
        # Use AiZynthFinder for retrosynthesis
        try:
            from syndirella.aizynth.AiZynthManager import AiZynthManager
            aizynth_search = AiZynthManager()
            routes: List[Dict[str, List[Dict[str, str]]]] = aizynth_search.perform_route_search(
                target_smiles=scaffold,
                matching_strategy='best_overall'
            )
            # Convert AiZynthFinder format to expected format
            routes = _convert_aizynth_routes(routes)
        except ImportError:
            logger.error("AiZynthFinder not available. Falling back to Manifold.")
            postera = Postera()
            routes: List[Dict[str, List[Dict[str, str]]]] | None = postera.perform_route_search(scaffold)
            if routes is None:
                logger.critical(f"API retrosynthesis query failed for {scaffold}.")
                raise APIQueryError(message=f"API retrosynthesis query failed for {scaffold}.", smiles=scaffold)
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


def process_df(df: pd.DataFrame, retro_tool: RetrosynthesisTool = DEFAULT_RETROSYNTHESIS_TOOL):
    """
    Process the input DataFrame and create output df with retrosynthesis information.
    
    Args:
        df: Input DataFrame with 'smiles' column
        retro_tool: Retrosynthesis tool to use
    """
    logger.info(f"Processing DataFrame with len {len(df)} using {retro_tool}")
    route_infos = []
    for i, scaffold in enumerate(df['smiles']):
        logger.info(f"Processing scaffold {i+1}/{len(df)}: {scaffold}")
        route_info: pd.DataFrame = retro_search(scaffold, retro_tool)
        route_infos.append(route_info)
    # format the DataFrame
    route_info_df = pd.concat(route_infos)
    merged_df = df.merge(route_info_df, on='smiles')
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

    df: pd.DataFrame = process_df(df, retro_tool)

    saved_path: str = save_df(df, output_dir, csv_path)

    logger.info(f"Saved DataFrame to {saved_path}")

    logger.info('Justretroquery execution completed successfully.')
