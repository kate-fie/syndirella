#!/usr/bin/env python3
"""
CLI script to run the Syndirella pipeline or justretroquery 🏁
"""
import argparse
import cProfile
import importlib.util
import logging
import os
import sys
import traceback
from typing import Dict, Any


def setup_logging(level=logging.INFO):
    """
    Set up logging configuration
    """
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )


def load_pipeline_module(syndirella_base_path: str, module_relative_path: str):
    """
    Dynamically load the pipeline module from a specified path.
    """
    pipeline_module_path = os.path.join(syndirella_base_path, module_relative_path)
    spec = importlib.util.spec_from_file_location("syndirella.pipeline", pipeline_module_path)
    pipeline = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pipeline)
    return pipeline


def config_parser(syndirella_base_path: str):
    """
    Configure command-line argument parsing.
    """
    parser = argparse.ArgumentParser(prog="syndirella",
                                     description="Run the Syndirella pipeline with specified configurations.",
                                     epilog=f"Syndirella is installed at {syndirella_base_path} \n",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str, required=True, help="Input .csv file path for the pipeline.")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output directory for the pipeline results.")
    parser.add_argument('-t', '--templates', type=str, required=False,
                        help="Absolute path to a directory containing the "
                             "template(s).")
    parser.add_argument('--hits_path', type=str, required=False,
                        help="Absolute path to hits_path for placements (.sdf or .mol).")
    parser.add_argument('--metadata', type=str, required=False, help="Absolute path to metadata for placements.")
    parser.add_argument('--products', type=str, required=False, help="Absolute path to products for placements.")
    parser.add_argument('--batch_num', type=int, default=10000, help="Batch number for processing.")
    parser.add_argument('--manual', action='store_true', help="Use manual routes for processing.")
    parser.add_argument('--only_scaffold_place', action='store_true',
                        help="Only place scaffolds. Do not continue to elaborate.")
    parser.add_argument('--scaffold_place_num', type=int, default=5,
                        help="Number of times to attempt scaffold placement.")
    parser.add_argument('--profile', action='store_true', help="Run the pipeline with profiling.")
    parser.add_argument('--atom_diff_min', type=int, default=0,
                        help="Minimum atom difference between elaborations and scaffold to keep.")
    parser.add_argument('--atom_diff_max', type=int, default=10,
                        help="Maximum atom difference between elaborations and scaffold to keep.")
    parser.add_argument('--long_code_column', type=str, default='Long code',
                        help="Column name for long code in metadata csv to match to SDF name."
                             " The column can contain a substring for the SDF name. ")
    parser.add_argument('--just_retro', action='store_true', help="Only run retrosynthesis querying of scaffolds.")
    parser.add_argument('--no_scaffold_place', action='store_true',
                        help="Do not place scaffolds initially before elaborating, "
                             "immediately start elaboration process.")
    parser.add_argument('--elab_single_reactant', action='store_true',
                        help="Only elaborate one reactant per elaboration series. "
                             "Warning: Functionality only provided for single step reactions.")
    return parser


def run_pipeline(settings: Dict[str, Any], pipeline):
    """
    Run the pipeline with the given settings.
    """
    logging.info('Running the pipeline...')

    pipeline.run_pipeline(settings)

    logging.info('Pipeline execution completed successfully.')


def run_justretroquery(settings: Dict[str, Any], justretroquery):
    """
    Run the justretroquery pipeline with the given settings.
    """
    logging.info('Running justretroquery...')

    justretroquery.run_justretroquery(settings)

    logging.info('Justretroquery execution completed successfully.')


def main():
    """
    Main entry point for the CLI.
    """
    setup_logging()
    syndirella_base_path = os.path.dirname(importlib.util.find_spec('syndirella').origin)

    parser = config_parser(syndirella_base_path)
    args = parser.parse_args()

    # check for MANIFOLD API key
    if not os.environ.get('MANIFOLD_API_KEY'):
        logging.error("MANIFOLD_API_KEY environment variable not set.")
        sys.exit(1)

    # check for MANIFOLD API key
    if not os.environ.get('MANIFOLD_API_URL'):
        logging.error("MANIFOLD_API_URL environment variable not set.")
        sys.exit(1)

    # Convert argparse Namespace to dictionary
    settings = vars(args)

    if settings['just_retro']:
        # Load the justretroquery module
        justretroquery = load_pipeline_module(syndirella_base_path, 'justretroquery.py')
        # Run the pipeline
        try:
            run_justretroquery(settings, justretroquery)
            sys.exit(0)
        except Exception as e:
            tb = traceback.format_exc()
            logging.error(f"An error occurred during justretroquery execution: {tb}")
            sys.exit(1)

    # Load the pipeline module
    pipeline = load_pipeline_module(syndirella_base_path, 'pipeline.py')

    if settings['profile']:
        try:
            profiler = cProfile.Profile()
            profiler.enable()
            run_pipeline(settings, pipeline)
            profiler.disable()
            profiler.print_stats(sort='time')
            print('\n\n')
            profiler.print_stats(sort='cumtime')
            sys.exit(0)
        except Exception as e:
            tb = traceback.format_exc()
            logging.error(f"An error occurred during pipeline execution: {tb}")
            sys.exit(1)

    # Run the pipeline
    try:
        run_pipeline(settings, pipeline)
    except Exception as e:
        tb = traceback.format_exc()
        logging.error(f"An error occurred during pipeline execution: {tb}")
        sys.exit(1)


if __name__ == '__main__':
    main()
