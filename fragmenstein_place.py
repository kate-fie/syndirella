import ast

from rdkit.ML.Cluster import Butina
from rdkit.Chem import rdMolDescriptors as rdmd
from rdkit.Chem import Descriptors

from typing import List, Dict, Any, Optional
import operator, os, re, logging, random, time, argparse, string, itertools, json, contextlib, requests
from warnings import warn
from pathlib import Path

import pyrosetta
import pyrosetta_help as ph

import pandas as pd
import pandera.typing as pdt
from pandarallel import pandarallel


from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import AllChem, Draw, PandasTools, BRICS
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdFingerprintGenerator as rdfpg
from rdkit.Chem.rdfiltercatalog import FilterCatalogParams, FilterCatalog, FilterCatalogEntry


Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

from fragmenstein import Igor, Victor, Laboratory, Monster
from fragmenstein.laboratory.validator import place_input_validator

pandarallel.initialize()

SUPRESSED_EXCEPTION = Exception  # debug purposes...
DEFAULT_WEIGHTS = {"N_rotatable_bonds": 3, "\u2206\u2206G": 3, "interaction_uniqueness_metric": -20, "N_unconstrained_atoms": 0.5, "N_constrained_atoms": -0.2, "N_interactions": -5, "N_interactions_lost": 10, "max_hit_Tanimoto": -1, "N_PAINS": 20, "strain_per_HA": 1}
# ------------------------------------------------------

def set_up(output, cutoff, quick, suffix, **kwargs):
    os.makedirs(output, exist_ok=True)
    Victor.work_path = output
    Laboratory.work_path = output
    Victor.monster_throw_on_discard = True  # stop this merger if a fragment cannot be used.
    Victor.monster_joining_cutoff = cutoff  # Å
    Victor.quick_reanimation = quick  # for the impatient
    Victor.error_to_catch = Exception  # stop the whole laboratory otherwise
    Victor.enable_stdout(logging.CRITICAL)
    Victor.enable_logfile(os.path.join(output, f'{suffix}.log'), logging.ERROR)


def config_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--template', help='Template PDB file', required=True)
    parser.add_argument('-i', '--hits', help='Hits file (.mol or .sdf)', required=True)
    parser.add_argument('-a', '--analog_path', help='Analogs csv file', required=True)
    parser.add_argument('-o', '--output', help='Output folder', default='output')
    parser.add_argument('-r', '--ranking', help='Ranking method', default='∆∆G')
    parser.add_argument('-c', '--cutoff', help='Joining cutoff', default=5)
    parser.add_argument('-q', '--quick', help='Quick reanimation', default=False)
    parser.add_argument('-s', '--suffix', help='Suffix for output files', default='')
    parser.add_argument('-n', '--n_cores', help='Number of cores', default=8, type=int)
    parser.add_argument('-m', '--combination_size', help='Number of hits to combine in one step', default=2, type=int)
    parser.add_argument('-k', '--top_mergers', help='Max number of mergers to followup up on', default=500, type=int)
    parser.add_argument('-e', '--timeout', help='Timeout for each merger', default=240, type=int)
    parser.add_argument('-x', '--max_tasks', help='Max number of combinations to try', default=0, type=int)
    parser.add_argument('-z', '--blacklist', help='Blacklist file', default='')
    parser.add_argument('-j', '--weights', help='JSON weights file', default='')
    return parser

# ------------------------------------------------------


def replace_hits(pdbblock, hits, n_cores, timeout, suffix, **settings):
    # place themselves for a ∆∆G score
    lab = Laboratory(pdbblock=pdbblock, covalent_resi=None, run_plip=True)
    selfies = pd.DataFrame([dict(name=hit.GetProp('_Name'),
                                 hits=[hit],
                                 smiles=Chem.MolToSmiles(hit)
                                 ) for hit in hits])
    replacements: pd.DataFrame = lab.place(place_input_validator(selfies), n_cores=n_cores, timeout=timeout)
    fix_intxns(replacements)
    replacements['bleached_name'] = replacements['name']
    replacements['name'] = replacements.hit_mols.apply(lambda ms: ms[0].GetProp('_Name'))
    replacements.to_pickle(f'fragmenstein_hit_replacements{suffix}.pkl.gz')
    # replacements.to_csv(f'fragmenstein_hit_replacements{suffix}.csv')
    return replacements

def UFF_Gibbs(mol):
    # free energy cost of bound conformer
    if not isinstance(mol, Chem.Mol) or mol.GetNumHeavyAtoms() == 0:
        return float('nan')
    with contextlib.suppress(SUPRESSED_EXCEPTION):
        AllChem.SanitizeMol(mol)
        # this is actually UFF
        copy = Chem.Mol(mol)
        return Monster.MMFF_score(None, mol, True)
    return float('nan')

# ------------------------------------------------------

def place(similars, pdbblock, n_cores, timeout, suffix, **settings) -> pd.DataFrame:
    """

    :param similars: Needs to have column, 'name', 'hits', 'smiles'
    :param pdbblock:
    :param n_cores:
    :param timeout:
    :param suffix:
    :param settings:
    :return:
    """
    lab = Laboratory(pdbblock=pdbblock, covalent_resi=None, run_plip=True)
    lab.work_path = settings['output']
    placements: pd.DataFrame = lab.place(place_input_validator(similars), n_cores=n_cores, timeout=timeout)
    placements.loc[(placements['∆∆G'] > -1) & (placements.outcome == 'acceptable'), 'outcome'] = 'weak'
    placements['unminimized_mol'] = placements.unminimized_mol.fillna(Chem.Mol())
    fix_intxns(placements)
    placements.to_pickle(f'fragmenstein_placed{suffix}.pkl.gz')
    placements.to_csv(f'fragmenstein_placed{suffix}.csv')
    print(placements.outcome.value_counts())
    return placements

def fix_intxns(df):
    intxn_names = [c for c in df.columns if isinstance(c, tuple)]
    for intxn_name in intxn_names:
        df[intxn_name] = df[intxn_name].fillna(0).astype(int)

def write(placements, suffix):
    valids = placements.loc[placements.outcome == 'acceptable'].sort_values('∆∆G')
    valids['combined_name'] = valids['name']
    valids['name'] = valids['enamine_name']
    valids.to_pickle(f'fragmenstein_acceptables{suffix}.pkl.gz')
    valids.path.apply(Path.exists).value_counts()

def correct_weaklings(hit_replacements, df):
    df['hit_names'] = df.hit_mols \
        .apply(lambda v: v if isinstance(v, list) else []) \
        .apply(lambda v: [m.GetProp('_Name') for m in v])
    dG_mapping = hit_replacements.set_index('name')['∆∆G'].to_dict()
    dG_mapping.update({k.replace('-', '_'): v for k, v in dG_mapping.items()})
    get_lowest = lambda names: min([0] + [dG_mapping.get(name, 0) for name in names])
    df['lowest_hit_∆∆G'] = df.hit_names.apply(get_lowest)
    worseness_mask = (df['∆∆G'] > df['lowest_hit_∆∆G'] * 0.8) & (
            df.outcome == 'acceptable')
    df.loc[worseness_mask, 'outcome'] = 'weaker'

def add_hits(analogs, **settings):
    analogs['hits'] = analogs['smiles'].apply(lambda x: settings['hits'])
    return analogs


def core_ops(analogs, **settings):
    # add hits column
    analogs = add_hits(analogs, **settings)
    placements: pd.DataFrame = place(analogs, **settings)
    return placements

# ------------------------------------------------------

if __name__ == '__main__':
    parser = config_parser()
    # load
    settings: Dict[str, Any] = vars(parser.parse_args())
    set_up(**settings)
    # load hits either from mol or sdf
    if os.path.splitext(settings['hits'])[1] == '.mol':
        print('This is a mol file')
        hits: List[Chem.Mol] = [Chem.MolFromMolFile(settings['hits'].strip())]
    else:
        with Chem.SDMolSupplier(settings['hits'].strip()) as sd:
            # hitdex: Dict[str, Chem.Mol] = {mol.GetProp('_Name'): mol for mol in sd}
            hits: List[Chem.Mol] = list(sd)
    settings['hits'] = hits
    if settings['blacklist']:
        with open(settings['blacklist'].strip()) as fh:
            settings['blacklist'] = [line.strip() for line in fh.readlines()]
    else:
        settings['blacklist'] = []
    print(f'N hits: {len(hits)}')
    with open(settings['template'].strip()) as fh:
        pdbblock = fh.read()
    settings['pdbblock'] = pdbblock
    if settings['weights']:
        with open(settings['weights'].strip()) as fh:
            settings['weights'] = json.load(fh)
    else:
        settings['weights'] = DEFAULT_WEIGHTS
    # self
    hit_replacements: pd.DataFrame = replace_hits(**settings)
    # run
    max_tasks = settings['max_tasks']
    hitnames = [h.GetProp('_Name') for h in hits]
    all_names = list(map('-'.join, itertools.permutations(hitnames, settings['combination_size'])))
    analogs = pd.read_csv(settings['analog_path'])
    if max_tasks == 0 or max_tasks > len(all_names):
        core_ops(analogs, **settings)
        exit()
    base_suffix = settings['suffix']
    all_placements = pd.DataFrame()
    letters = iter(string.ascii_uppercase)
    for i in range(0, len(all_names) + max_tasks, max_tasks):
        settings['suffix'] = base_suffix + next(letters)
        placements: pd.DataFrame = core_ops(analogs, **settings)
        settings['blacklist'] += all_names[i:i + max_tasks]
        all_placements = pd.concat([all_placements, placements], ignore_index=True)
    correct_weaklings(hit_replacements, all_placements)
    settings['suffix'] = base_suffix
    all_placements.to_pickle(f'fragmenstein_placed{base_suffix}.pkl.gz')
    all_placements.to_pickle(f'fragmenstein_placed{base_suffix}.pkl.gz')
    all_placements.to_csv(f'fragmenstein_placed{base_suffix}.csv')
