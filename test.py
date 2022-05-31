#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 13:29:21 2021

@author: ruben
"""
from chemUtils.substructure import splitExtendedMolBySubstruct
from rdkit import Chem
from rdkit.Chem import AllChem


data_examples_from_bbs=[
  dict(ori_smi="CC1OCCC1C(=O)NCCc1ccccc1", toAdd_smi="Cn1cccc1CNC(=O)C1CCC(N)CC1",
       ori_attachmentIdx=13, frag_attachmentIdx=31, post_attachmentIdx=15),
  dict(ori_smi="CC(C1=CC=CC=C1)N2CC(CC2=O)C(=O)O", toAdd_smi="Fc1ccc(CNCc2ccc3c(c2)OCO3)cc1",
       ori_attachmentIdx=16, frag_attachmentIdx=23, post_attachmentIdx=None),
]

data_examples_from_fragment_full=[
  dict(ori_smi="COC1C=C(CO)C=CC=1O", full_smi="Clc1cccc(F)c1CNCc1ccc(O)c(OC)c1"),
  dict(ori_smi="C1=CC(=C(C(=C1)Cl)CN)F", full_smi="Clc1cccc(F)c1CNCc1ccc(O)c(OC)c1"),

  dict(ori_smi="C1(C)C=C(C#N)C=NC=1", full_smi="Cc1ccc(N(c2cncc(C#N)c2)C(N)C2CC2)cc1"),
  dict(ori_smi="C(CCCC1C=CC=CC=1)=O", full_smi="Nc1ccc(COC(=O)CCCc2ccccc2)cc1"),

]

def _getDataExample_fromComponents(ori_smi, toAdd_smi, ori_attachmentIdx, frag_attachmentIdx, post_attachmentIdx ):

  ori_mol = Chem.MolFromSmiles(ori_smi)
  toAdd_mol = Chem.MolFromSmiles(toAdd_smi)
  combined_mol = Chem.CombineMols(ori_mol, toAdd_mol)
  combined_mol = Chem.RWMol(combined_mol)
  # from molPlot import plotMols
  # plotMols([ori_mol, toAdd_mol, combined_mol], show3D=False)
  combined_mol.AddBond(ori_attachmentIdx,frag_attachmentIdx,Chem.BondType.SINGLE)
  combined_mol = combined_mol.GetMol()
  combined_mol = AllChem.AssignBondOrdersFromTemplate(combined_mol, Chem.MolFromSmiles(Chem.MolToSmiles(combined_mol)))
  AllChem.EmbedMolecule(combined_mol)

  (oriMol, ori_attachmentIdx), _, extendedMol_attachmentIdx = splitExtendedMolBySubstruct(combined_mol, ori_mol)
  # from molPlot import plotMols
  # plotMols([ori_mol, toAdd_mol, combined_mol], show3D=False)
  return post_attachmentIdx, oriMol, combined_mol

def _getDataExample_fromFragmentAndFinal(ori_smi, full_smi ):

  ori_mol = Chem.MolFromSmiles(ori_smi)
  combined_mol = Chem.MolFromSmiles(full_smi)

  AllChem.EmbedMolecule(combined_mol)

  (oriMol, ori_attachmentIdx), _, extendedMol_attachmentIdx = splitExtendedMolBySubstruct(combined_mol, ori_mol)
  # from molPlot import plotMols
  # plotMols([ori_mol, toAdd_mol, combined_mol], show3D=False)

  return extendedMol_attachmentIdx, oriMol, combined_mol

def _test_searchDirectAnalogues():
  from wholeMoleculePipeline import searchDirectAnalogues
  attachmentIdxs, ori_mol, modified_mol = _getDataExample_fromComponents(**data_examples_from_bbs[0])
  results = searchDirectAnalogues(modified_mol)
  # print(results)
  for score, smi, similarity in results:
    print(smi, score, similarity)
  print( Chem.MolToSmiles(modified_mol), "TARGET")

def _test_searchPicewiseAnalogues1():
  from wholeMoleculePipeline import searchPicewiseAnalogues
  ori_attachmentIdx, ori_mol, modified_mol = _getDataExample_fromComponents(**data_examples_from_bbs[0])
  results = searchPicewiseAnalogues(ori_mol, modified_mol, ori_attachmentIdx)
  print(results)

def _test_searchPicewiseAnalogues2():
  from wholeMoleculePipeline import searchPicewiseAnalogues
  ori_attachmentIdx, ori_mol, modified_mol = _getDataExample_fromFragmentAndFinal(**data_examples_from_fragment_full[2])
  results = searchPicewiseAnalogues(ori_mol, modified_mol)
  print(results)

def _test_figure_out_attachment(): #TODO
  raise NotImplementedError()



# test_searchDirectAnalogues()
# _test_searchPicewiseAnalogues1()
_test_searchPicewiseAnalogues2()


