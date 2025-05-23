#!/usr/bin/env python3
"""
syndirella.run_pipeline.py

This script contains the main pipeline for syndirella.
"""
import datetime
import logging
import time
import traceback
from typing import List, Tuple, Dict

import pandas as pd
from rdkit.Chem import inchi
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

import syndirella.check_inputs as check_inputs
import syndirella.fairy as fairy
import syndirella.structure_outputs as structure_outputs
from syndirella.Cobbler import Cobbler
from syndirella.error import *
from syndirella.route.CobblersWorkshop import CobblersWorkshop
from syndirella.slipper.Slipper import Slipper
from syndirella.slipper.SlipperFitter import SlipperFitter

logger = logging.getLogger(__name__)


def assert_scaffold_placement(scaffold: str,
                              template_path: str,
                              hits_path: str,
                              hits_names: List[str],
                              output_dir: str,
                              scaffold_place_num: int
                              ) -> Dict[Chem.Mol, str]:
    """
    Assert that the scaffold can be placed for any stereoisomers. If not, raise an error.
    """
    scaffold_mol = Chem.MolFromSmiles(scaffold)
    if not scaffold_mol:
        logger.critical(f"Could not create a molecule from the smiles {scaffold}.")
        raise MolError(smiles=scaffold)
    # enumerate stereoisomers
    opts = StereoEnumerationOptions(unique=True)
    isomers = list(EnumerateStereoisomers(scaffold_mol, options=opts))
    slipper_fitter = SlipperFitter(template_path, hits_path, hits_names, output_dir)
    placements: Dict[Chem.Mol, str | None] = {}
    for i, isomer in enumerate(isomers):
        scaffold_name: str = f'scaffold-{chr(65 + i)}'
        can_be_placed: str | None = slipper_fitter.check_scaffold(scaffold=isomer,
                                                                  scaffold_name=scaffold_name,
                                                                  scaffold_place_num=scaffold_place_num)  # path to scaffold if successful
        placements[isomer] = can_be_placed  # absolute path to minimised.mol scaffold, checked to exist
    if not any(placements.values()):
        logger.critical(f"Scaffold {scaffold} could not be placed successfully.")
        raise ScaffoldPlacementError(smiles=scaffold)
    return placements


def elaborate_from_cobbler_workshops(cobbler_workshops: List[CobblersWorkshop],
                                     template_path: str,
                                     hits_path: str,
                                     hits: List[str],
                                     batch_num: int,
                                     csv_path: str,
                                     output_dir: str,
                                     scaffold_placements: Dict[Chem.Mol, str],
                                     additional_info=None):
    """
    Does elaboration once the cobbler workshops are created.
    """
    if additional_info is None:
        additional_info = []
    for workshop in cobbler_workshops:
        try:
            slipper = None
            final_library = workshop.get_final_library()
            if final_library is None:
                logger.warning(f"Could not get the final library for compound {workshop.product}. Skipping...")
                continue
            slipper = Slipper(library=final_library, template=template_path, hits_path=hits_path, hits_names=hits,
                              batch_num=batch_num, atoms_ids_expansion=None, additional_info=additional_info,
                              scaffold_placements=scaffold_placements)
            slipper.get_products()
            slipper.place_products()
            slipper.write_products_to_hippo()  # only write at the end after placement, to get correct route_uuid
            slipper.clean_up_placements()
        except Exception as e:
            tb = traceback.format_exc()
            logger.critical(f"Error elaborating compound {workshop.product}. {tb}")
            structure_outputs.structure_pipeline_outputs(error=e,
                                                         csv_path=csv_path,
                                                         output_dir=output_dir,
                                                         workshop=workshop if workshop is not None else None,
                                                         slipper=slipper if slipper is not None else None)
            continue
        structure_outputs.structure_pipeline_outputs(csv_path=csv_path,
                                                     output_dir=output_dir,
                                                     error=None,
                                                     workshop=workshop,
                                                     slipper=slipper)


def start_elaboration(product: str,
                      template_path: str,
                      hits_path: str,
                      hits: List[str],
                      output_dir: str,
                      scaffold_place_num: int,
                      scaffold_place: bool) -> Tuple[float, Dict[Chem.Mol, str | None] | Dict]:
    """
    Starts the elaboration of a single compound.
    """
    logger.info(f'Starting: {product} | {inchi.MolToInchiKey(Chem.MolFromSmiles(product))}')
    start_time = time.time()
    if not scaffold_place:
        logger.info(f"Skipping initial scaffold placement...")
        scaffold_placements: Dict[Chem.Mol, str | None] = {}
    else:
        logger.info(f"Placing scaffold...")
        scaffold_placements: Dict[Chem.Mol, str | None] = assert_scaffold_placement(scaffold=product,
                                                                                    template_path=template_path,
                                                                                    hits_path=hits_path,
                                                                                    hits_names=hits,
                                                                                    output_dir=output_dir,
                                                                                    scaffold_place_num=scaffold_place_num)
    return start_time, scaffold_placements


def elaborate_compound_with_manual_routes(product: str,
                                          reactants: List[Tuple[str, str]],
                                          reaction_names: List[str],
                                          num_steps: int,
                                          hits: List[str],
                                          template_path: str,
                                          hits_path: str,
                                          batch_num: int,
                                          output_dir: str,
                                          csv_path: str,
                                          atom_diff_min: int,
                                          atom_diff_max: int,
                                          scaffold_place_num: int,
                                          only_scaffold_place: bool,
                                          scaffold_place: bool,
                                          elab_single_reactant: bool,
                                          additional_info=None):
    """
    This function is used to elaborate a single compound using a manually defined route.
    """
    start_time, scaffold_placements = start_elaboration(product=product, template_path=template_path,
                                                        hits_path=hits_path,
                                                        hits=hits, output_dir=output_dir,
                                                        scaffold_place_num=scaffold_place_num,
                                                        scaffold_place=scaffold_place)
    if not only_scaffold_place:  # continue elaboration
        workshop = CobblersWorkshop(product=product, reactants=reactants, reaction_names=reaction_names,
                                    num_steps=num_steps, output_dir=output_dir, filter=False,
                                    id=fairy.generate_inchi_ID(product, isomeric=False), atom_diff_min=atom_diff_min,
                                    atom_diff_max=atom_diff_max, elab_single_reactant=elab_single_reactant)
        cobblers_workshops = [workshop]
        alternative_routes: List[CobblersWorkshop] | None = workshop.get_additional_routes(edit_route=True)
        if alternative_routes is not None:
            cobblers_workshops = [workshop] + alternative_routes
        elaborate_from_cobbler_workshops(cobbler_workshops=cobblers_workshops, template_path=template_path,
                                         hits_path=hits_path, hits=hits, batch_num=batch_num,
                                         additional_info=additional_info, csv_path=csv_path, output_dir=output_dir,
                                         scaffold_placements=scaffold_placements)
        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.info(
            f"Finished Syndirella 👑 pipeline for compound {product} | {inchi.MolToInchiKey(Chem.MolFromSmiles(product))} "
            f"after {datetime.timedelta(seconds=elapsed_time)}")
        logger.info("")


def elaborate_compound_full_auto(product: str,
                                 hits: List[str],
                                 template_path: str,
                                 hits_path: str,
                                 batch_num: int,
                                 output_dir: str,
                                 csv_path: str,
                                 atom_diff_min: int,
                                 atom_diff_max: int,
                                 scaffold_place_num: int,
                                 only_scaffold_place: bool,
                                 scaffold_place: bool,
                                 elab_single_reactant: bool,
                                 additional_info=None):
    """
    This function is used to elaborate a single compound.
    """
    start_time, scaffold_placements = start_elaboration(product=product, template_path=template_path,
                                                        hits_path=hits_path,
                                                        hits=hits, output_dir=output_dir,
                                                        scaffold_place_num=scaffold_place_num,
                                                        scaffold_place=scaffold_place)
    if not only_scaffold_place:  # continue elaboration
        cobbler = Cobbler(scaffold_compound=product,
                          output_dir=output_dir, atom_diff_min=atom_diff_min,
                          atom_diff_max=atom_diff_max,
                          elab_single_reactant=elab_single_reactant)  # check that output_dirs can be made for different routes
        cobbler_workshops: List[CobblersWorkshop] = cobbler.get_routes()
        elaborate_from_cobbler_workshops(cobbler_workshops=cobbler_workshops, template_path=template_path,
                                         hits_path=hits_path, hits=hits, batch_num=batch_num,
                                         additional_info=additional_info, csv_path=csv_path, output_dir=output_dir,
                                         scaffold_placements=scaffold_placements)
        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.info(
            f"Finished Syndirella 👑 pipeline for compound {product} | {inchi.MolToInchiKey(Chem.MolFromSmiles(product))} "
            f"after {datetime.timedelta(seconds=elapsed_time)}")
        logger.info("")


#######################################

def run_pipeline(settings: Dict):
    """
    Run the whole syndirella pipeline! 👑
    """

    def process_row(row: pd.Series, manual_routes: bool, only_scaffold_place: bool, scaffold_place: bool):
        additional_info: dict = check_inputs.format_additional_info(row, additional_columns)
        template_path: str = check_inputs.get_template_path(
            template_dir=template_dir,
            template=row['template'],
            metadata_path=metadata_path
        )

        hits: List[str] = check_inputs.get_exact_hit_names(row=row, metadata_path=metadata_path,
                                                           hits_path=hits_path)

        # logger.warning("Assuming hit names in SDF exactly match input CSV")
        # hits = [v for k,v in row.to_dict().items() if k.startswith("hit")]

        try:
            if manual_routes:
                reactants, reaction_names, num_steps = check_inputs.format_manual_route(row)
                elaborate_compound_with_manual_routes(
                    product=row['smiles'],
                    reactants=reactants,
                    reaction_names=reaction_names,
                    num_steps=num_steps,
                    hits=hits,
                    template_path=template_path,
                    hits_path=hits_path,
                    batch_num=batch_num,
                    output_dir=output_dir,
                    additional_info=additional_info,
                    csv_path=csv_path,
                    atom_diff_min=atom_diff_min,
                    atom_diff_max=atom_diff_max,
                    scaffold_place_num=scaffold_place_num,
                    only_scaffold_place=only_scaffold_place,
                    scaffold_place=scaffold_place,
                    elab_single_reactant=elab_single_reactant
                )
            else:
                elaborate_compound_full_auto(
                    product=row['smiles'],
                    hits=hits,
                    template_path=template_path,
                    hits_path=hits_path,
                    batch_num=batch_num,
                    output_dir=output_dir,
                    additional_info=additional_info,
                    csv_path=csv_path,
                    atom_diff_min=atom_diff_min,
                    atom_diff_max=atom_diff_max,
                    scaffold_place_num=scaffold_place_num,
                    only_scaffold_place=only_scaffold_place,
                    scaffold_place=scaffold_place,
                    elab_single_reactant=elab_single_reactant
                )
        except Exception as e:
            tb = traceback.format_exc()
            logger.critical(f"Error elaborating compound {row['smiles']}. {tb}")
            structure_outputs.structure_pipeline_outputs(
                error=e,
                csv_path=csv_path,
                output_dir=output_dir,
                smiles=row['smiles'],
                template_path=template_path,
                hits=hits,
                additional_info=additional_info
            )

    # set required variables
    try:
        additional_columns: List[str] = ['compound_set']
        metadata_path: str = settings['metadata']
        template_dir: str = settings['templates']
        hits_path: str = settings['hits_path']
        output_dir: str = settings['output']
        batch_num: int = settings['batch_num']
        csv_path: str = settings['input']
        atom_diff_min: int = settings['atom_diff_min']
        atom_diff_max: int = settings['atom_diff_max']
        scaffold_place_num: int = settings['scaffold_place_num']
        long_code_column: str = settings['long_code_column']
    except KeyError as e:
        logger.critical(f"Missing critical argument to run pipeline: {e}")

    # set optional variables
    try:
        manual_routes: bool = settings['manual']
    except KeyError:
        manual_routes = False

    try:
        only_scaffold_place: bool = settings['only_scaffold_place']
        # If only_scaffold_place is True, only place scaffolds and do not continue to elaborate
        if only_scaffold_place:
            logger.info(f"Only placing scaffolds!")
    except KeyError:
        only_scaffold_place = False
        # Log pipeline type

    try:
        # if no_scaffold_place is True, do not place scaffolds
        scaffold_place: bool = not settings['no_scaffold_place']
        if not scaffold_place:
            logger.warning(f"Skipping initial scaffold placement! Immediately starting elaboration process.")
    except KeyError:
        scaffold_place = True

    if not only_scaffold_place:
        logger.info(f"Running the pipeline with {'manual' if manual_routes else 'full auto'} routes.")

    try:
        elab_single_reactant: bool = settings['elab_single_reactant']
        if elab_single_reactant:
            logger.info(f"'--elab_single_reactant' set. Only elaborating a single reactant per elaboration series.")
    except KeyError:
        elab_single_reactant = False

    # Validate inputs
    check_inputs.check_pipeline_inputs(
        csv_path=csv_path,
        template_dir=template_dir,
        hits_path=hits_path,
        metadata_path=metadata_path,
        additional_columns=additional_columns,
        manual_routes=manual_routes,
        long_code_column=long_code_column
    )

    # Load data
    df = pd.read_csv(csv_path)

    # Process each row in the DataFrame
    for index, row in df.iterrows():
        process_row(row=row, manual_routes=manual_routes, only_scaffold_place=only_scaffold_place,
                    scaffold_place=scaffold_place)

    logger.info("Pipeline complete.")
