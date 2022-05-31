#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:56:02 2021

@author: ruben
"""
import sys
from typing import TextIO

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

import config
from config import config

if config.SUCOS_DIR not in sys.path:
  sys.path.append(config.SUCOS_DIR)

from sucos import SuCOS


def calc_energy(mol: Chem.Mol) -> float:
  """
  Function to calculate the energy of the conformer 0 of a molecule
  :param mol: embedded molecule
  :return: energy of the molecule
  :rtype: float
  """
  try:
    mol_energy = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
  except RuntimeError:
    mol_energy = np.inf
  return mol_energy

def calc_unconstrained_energy(og_mol: Chem.Mol, n_conformers=config.N_CONFORMERS_FOR_CONSTRAIN_EMBEDDING_FILTER) -> float:
  """
  Create ten unconstrained conformations for each molecule and calculate the energy.
  :param og_mol: the original molecule with given coordinates
  :return: the average of the unconstrained energies
  """
  og_mol = Chem.Mol(og_mol)
  og_mol.RemoveAllConformers()

  # try:
  confIds = list(AllChem.EmbedMultipleConfs(og_mol, numConfs=n_conformers, numThreads=config.MULTIEMBEDDING_N_THREADS))
  if len(confIds) == 0:
    return np.inf
  ok_energy_list = AllChem.UFFOptimizeMoleculeConfs(og_mol, numThreads=config.MULTIEMBEDDING_N_THREADS, maxIters=config.MULTIEMBEDDING_MAX_ITERS)
  unconstrained_energies = [ ok_energy[1] for ok_energy in ok_energy_list  ]
  # unconstrained_energies = [ ok_energy[1] for ok_energy in ok_energy_list if ok_energy[0]==0 ]

  if len(unconstrained_energies) == 0:
    avg = np.inf
  else:
    avg = sum(unconstrained_energies) / len(unconstrained_energies)
  return avg


def scorePairOfMols(originalMol, newMol):

  original_energy = calc_energy(originalMol) #TODO: handle excpetions
  new_energy = calc_energy(newMol)
  unconst_energy = calc_unconstrained_energy(newMol)

  score = 0
  epsilon = 1e-10
  energyFold_ori_new = abs(max(original_energy, new_energy)/(min(original_energy, new_energy))+epsilon)
  energyFold_ori_unconstraint = abs(max(original_energy, unconst_energy)/(min(original_energy, unconst_energy)+epsilon))

  if energyFold_ori_new > config.ORI_VS_NEW_ENERGY_TIMES_THR or \
     energyFold_ori_unconstraint > config.CONSTRAIN_VS_UNCONSTRAIN_ENERGY_TIMES_THR:
    score = np.nan #- np.inf
  score += SuCOS.computeShapeIOU(newMol, originalMol)
  score += SuCOS.computeFeatsIOU(newMol, originalMol)
  score *= 0.5
  return score, energyFold_ori_new, energyFold_ori_unconstraint


def scoreSDF(refMolFile:str, toScoreFile:str, outName:str = None):
  '''

  :param refMolFile:
  :param toScoreFile:
  :param outName:
  :return:
  '''
  from engine import myMap

  refMol = Chem.MolFromMolFile(refMolFile)
  suppl = Chem.SDMolSupplier(toScoreFile)

  def scoringFunction(newMol):
    return scorePairOfMols(refMol, newMol)

  scores_list = myMap(scoringFunction, suppl)

  if outName:
    raise NotImplementedError()
  else:
    for score in scores_list:
      score = list(score) + [ score[1]/score[2]]
      print(" ".join(map(lambda x: str(round(x, 3)),score)))

  return scores_list
if __name__ == "__main__":
    from argParseFromDoc import get_parser_from_function
    parser = get_parser_from_function(scoreSDF, args_optional=["outName"])
    args = parser.parse_args()
    scoreSDF(**vars(args))