"""
This file contains the ability to compare the constrained to average unconstrained energy of multiple conformations.

Author: Stephanie Wills
Created: 10-August-2022
"""
import argparse
import os
import sys
from typing import Tuple

from rdkit import Chem, RDLogger
from rdkit.Chem import Mol, AllChem
from config_filter import config_filter

RDLogger.DisableLog("rdApp.*")

def calc_energy(mol: Mol) -> float:
    """
    Calculate energy of molecule.
    """
    mol_energy = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
    return mol_energy


def calc_unconstrained_energy(og_mol: Mol, n_conf: int) -> float:
    """
    Calculate average energy of multiple unconstrained conformations of a molecule.
    """
    unconstrained_energies = []
    for i in range(n_conf):  # generate conformations and calculate energy
        mol = Chem.Mol(og_mol)
        AllChem.EmbedMolecule(mol, randomSeed=i*10)
        AllChem.UFFOptimizeMolecule(mol)
        e = calc_energy(mol)
        unconstrained_energies.append(e)
    # print(unconstrained_energies)
    # calculate the average of all the energies
    avg = sum(unconstrained_energies) / len(unconstrained_energies)
    # print(avg)
    return avg


class EnergyFilter:
    def __init__(self, mol):
        self.mol = mol

    def filter_smi(
        self,
        energy_threshold: float = config_filter.ENERGY_THRESHOLD,
        n_conf: int = config_filter.N_CONFORMATIONS,
    ) -> bool:
        try:
            const_energy = calc_energy(self.mol)
            print("constrained", const_energy)
            unconst_energy = calc_unconstrained_energy(self.mol, n_conf)
            print("unconstrained", unconst_energy)
            if (const_energy / unconst_energy) >= energy_threshold:
                return False
            else:
                return True
        except Exception as e:
            print(e)
            return False


def parse_args(args):
    parser = argparse.ArgumentParser(
        epilog="""
        python energy_filter.py --input_mol data/input.mol --energy_threshold 10 --n_conformations 50
        """
    )
    parser.add_argument(
        "-i",
        "--input_mol",
        required=True,
        help="input .mol file containing a molecule to filter",
    )
    parser.add_argument(
        "-e",
        "--energy_threshold",
        type=float,
        default=10.0,
        help="threshold for fold difference in energy of unconstrained versus constrained conformations",
    )
    parser.add_argument(
        "-n",
        "--n_conformations",
        type=int,
        default=50,
        help="number of unconstrained conformations to generate",
    )

    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    mol = Chem.MolFromMolFile(args.input_mol)
    if not mol:
        print(f"Failed to read molecule from {args.input_mol}")
        sys.exit(1)

    filter_obj = EnergyFilter(mol)
    result = filter_obj.filter_smi(args.energy_threshold, args.n_conformations)

    print(f"Molecule {'passes' if result else 'fails'} the energy filter.")


if __name__ == "__main__":
    main()
