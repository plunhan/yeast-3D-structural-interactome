import os
import pandas as pd 
import numpy as np 
from pathlib import Path 
from analysis_tools import correlation_ORC_GIPS, correlation_ORC_SemSim

def main():

    # threshold for distinguishing different "interfaces"
    threshold = 1

    # cut-off for distinguishing hub proteins
    hubCutOff = 2
    
    # reference interactome name
    interactome_name = 'Yeast'

    # parent directory of all data files
    dataDir = Path('../data')

    # parent directory of all external data files
    extDir = dataDir / 'external'

    # parent directory of all processed data files
    procDir = dataDir / 'processed'

    # parent directory of Aim 1 results (basic stats)
    aim1Dir = procDir / 'aim_1'

    # parent directory of Aim 2 results (correlation and confounding factors)
    aim2Dir = procDir / 'aim_2'

    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name

    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'

    # directory of semantic similarity files
    semSimDir = interactomeDir / 'gosim'

    # input data files
    referenceInteractome = procDir / 'yeast_reference_interactome.txt'
    structuralInteractome = modelBasedDir / 'structural_interactome.txt'
    degreeCountStrInt = aim1Dir / 'degree_count_structural_interactome.txt'
    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'
    structuralInteractomeWithInteractingResidues = aim2Dir / 'interacting_residues.pkl'
    geneticInteractionProfile_rel = procDir / 'sga_GI_rel_matrix.txt'
    geneticInteractionProfile_all = procDir / 'sga_GI_full_matrix.txt'
    semSimBP = semSimDir / 'fastsemsim_output_refPPIs_bp_SimGIC'
    semSimMF = semSimDir / 'fastsemsim_output_refPPIs_mf_SimGIC'
    semSimCC = semSimDir / 'fastsemsim_output_refPPIs_cc_SimGIC'

    # output data files

    corOverlapResCountGIPS = aim2Dir / 'overlapping_residue_count_vs_genetic_GIPS.txt'
    corOverlapResCountSemSim = aim2Dir / 'overlapping_residue_count_vs_semetic_similarity.txt'

    if not aim2Dir.exists():
        os.makedirs(aim2Dir)

    # correlation between overlapping residue count and genetic interaction profile similarity
    if not corOverlapResCountGIPS.is_file(): 
        correlation_ORC_GIPS(structuralInteractome, 
                             degreeCountStrInt, 
                             IDMappingFile, 
                             geneticInteractionProfile_rel, 
                             hubCutOff, 
                             corOverlapResCountGIPS)

    # correlation between overlapping residue count and GO-based semantic similarity
    if not corOverlapResCountSemSim.is_file():
        correlation_ORC_SemSim(structuralInteractome, 
                               degreeCountStrInt, 
                               IDMappingFile, 
                               geneticInteractionProfile_rel, 
                               semSimBP, 
                               semSimMF, 
                               semSimCC, 
                               hubCutOff, 
                               corOverlapResCountSemSim)

if __name__ == '__main__':
    main()