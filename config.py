#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:58:11 2021

@author: ruben
"""
import os

class Config(): #TODO: metaclass singleton
    ROOT_DIR = os.path.abspath(os.path.join(__file__, ".."))

    DUMMY_SYMBOL = "*"

    SUCOS_DIR = os.path.expanduser("/Users/kate_fieseler/PycharmProjects/chemUtils/chemUtils/scoring/_fragmentBasedScores/sucos.py")
    POSTERA_MANIFOLD_DIR = os.path.expanduser("/Users/kate_fieseler/PycharmProjects/postera")

    POSTERA_MANIFOLD_CACHE_FNAME = os.path.join(ROOT_DIR, "similarityCache/similarityCache.sqlite")
    POSTERA_MANIFOLD_SUPERSTRUCTURE_CACHE_FNAME = os.path.join(ROOT_DIR, "superstructureCache/superstructureCache.sqlite")

    POSTERA_MANIFOLD_VENDORS = [ "chemspace", "emolecules", "enamine_bb", "enamine_made", "enamine_real", "mcule", "mcule_ultimate", "molport", "wuxi_bb_screening", "wuxi_galaxi" ],
    POSTERA_MANIFOLD_PATENTS = []
    POSTERA_MANIFOLD_PATENTS_QUERY_THIRD_PARTIES = False

    SIMILARITY_SEARCH_THR = 0.3 #0.4
    STRUCTURAL_SCORE_THR = None #0.7
    POSITIONAL_MAPPING_MAX_DIST = 1.5 #Angstroms
    SMARTS_N_HOPS_FROM_ATTACHMENT = 2
    N_RETRY_QUERY=10

    N_CONFORMERS_FOR_CONSTRAIN_EMBEDDING_FILTER = 10
    CONSTRAIN_VS_UNCONSTRAIN_ENERGY_TIMES_THR = 10
    ORI_VS_NEW_ENERGY_TIMES_THR = 10


    MULTIEMBEDDING_N_THREADS = 5 # 0 means all
    MULTIEMBEDDING_MAX_ITERS = 300

    DATABASE_SEARCH_FUN = "databaseSearch.datatabaseSearch_postera.postera_similarity_search"

config = Config()