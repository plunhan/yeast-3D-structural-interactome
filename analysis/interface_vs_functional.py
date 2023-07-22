import os 
import pandas as pd 
import numpy as np 
from pathlib import Path 
from interface_tools import (degree_of_interface_overlap_vs_functional_similarity, 
                             binning_data, 
                             check_paralogs, 
                             generate_partner_evalue_distribution, 
                             degree_of_interface_overlap_vs_functional_similarity_new)
from confounding_factor_tools import (blast_interactome_protein_pairs, 
                                      generate_essentiality_summary, 
                                      check_essentiality, 
                                      check_Pfam_domains, 
                                      generate_structural_interactome_protein_sequences, 
                                      check_pfam_paralog, 
                                      check_pfam_at_interface, 
                                      produce_structural_similarity_whole_protein, 
                                      produce_structural_similarity_Pfam_domains, 
                                      produce_sequence_similarity_Pfam_domains, 
                                      produce_sequence_similarity_covered_interfaces, 
                                      produce_structural_similarity_extended_interfaces)

def main(): 

    # cut-off defining hub proteins
    hubCutOff = 2

    # cut-off for selecting target sequence alignment in blast
    alignment_cutoff = 1

    # cut-off for selecting valid TM-score in TM-align
    strSim_cutoff = 0.5
    
    # reference interactome name
    interactome_name = 'Yeast'

    # parent directory of all data files
    dataDir = Path('../data')

    # parent directory of all external data files
    extDir = dataDir / 'external'

    # parent directory of all processed data files
    procDir = dataDir / 'processed'

    # parent directory of Aim 1 results
    aim1Dir = procDir / 'aim_1'

    # parent directory of Aim 2 results
    aim2Dir = procDir / 'aim_2'

    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'

    # directory of semantic similarity files
    semSimDir = interactomeDir / 'gosim'

    # directory of PPI models
    ppiModelDir = modelBasedDir / 'ppi_models'

    # directory of Pfam scan files
    PfamDir = Path('../Pfam/PfamScan')

    # directory of TM-align files
    TMalignDir = Path('../TMalign/TMalign')

    # input files
    proteinSequence = extDir / 'UP000002311_559292.fasta'
    protSeqDict = procDir / 'yeast_reference_sequences.pkl'
    structuralInteractome = modelBasedDir / 'structural_interactome.txt'
    structuralInteractomeWithInteractingResidues = aim2Dir / 'interacting_residues.pkl'
    degreeStrInt = aim1Dir / 'degree_count_structural_interactome.txt'
    coexpression = extDir / '2010.Hughes00.flt.knn.avg.pcl'
    geneticInteractionProfile = procDir / 'GI_profile.txt'
    semanticSimilarity_mf = semSimDir / 'fastsemsim_output_refPPIs_mf_SimGIC'
    semanticSimilarity_bp = semSimDir / 'fastsemsim_output_refPPIs_bp_SimGIC'
    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'
    interactomeSequence = procDir / 'interactome_sequences.fasta'
    paralogs = extDir / 'Fares et al. 2013.xls'
    WGDs = extDir / 'Byrne and Wolfe 2005.xlsx'
    sga_GI_matrix = procDir / 'sga_GI_rel_matrix.txt'
    essentialGeneList = procDir / 'essential_genes.pkl'

    # output files
    degreeOverlapVsFunctionalSimilarity = aim2Dir / 'interface_overlap_vs_functional_similarity.txt'
    querySeqTempFile = aim2Dir / 'query_temp_seq.txt'
    subjectSeqTempFile = aim2Dir / 'subject_temp_seq.txt'
    blastTempFile = aim2Dir / 'blast_temp.txt'
    seqSimStrInt = aim2Dir / 'identities_positives_structural_interactome.txt'

    essentialitySummary = aim2Dir / 'essentiality_summary.txt'
    confounding_essentiality = aim2Dir / 'confounding_essentiality.txt'

    structuralInteractomeProteinSequence = aim2Dir / 'structural_interactome_protein_sequences.fasta'
    confounding_sequence_similarity_local = aim2Dir / 'confounding_sequence_similarity_local_case_study.txt'
    confounding_sequence_similarity_local_new = aim2Dir / 'confounding_sequence_similarity_local_new.txt'

    confounding_structural_similarity = aim2Dir / 'confounding_structural_similarity.txt'
    confounding_structural_similarity_local = aim2Dir / 'confounding_structural_similarity_local.txt'
    confounding_structural_similarity_local_new = aim2Dir / 'confounding_structural_similarity_local_new_50.txt'
    TempFile_1 = aim2Dir / 'tempfile1.txt'
    TempFile_2 = aim2Dir / 'tempfile2.txt'
    TMalign_file = aim2Dir / 'TMscore.txt'

    # Probe the effect of possible confounding factors

    # Essentiality
    if not essentialitySummary.is_file(): 
        generate_essentiality_summary(structuralInteractome, 
                                      degreeStrInt, 
                                      IDMappingFile, 
                                      essentialGeneList, 
                                      essentialitySummary)

    if not confounding_essentiality.is_file(): 
        check_essentiality(structuralInteractome, 
                           degreeStrInt, 
                           IDMappingFile, 
                           sga_GI_matrix, 
                           essentialGeneList, 
                           confounding_essentiality)

    # Sequence similarity
    # identities and positives of two genes
    if not structuralInteractomeProteinSequence.is_file(): 
        generate_structural_interactome_protein_sequences(proteinSequence, 
                                                          structuralInteractome, 
                                                          structuralInteractomeProteinSequence)

    if not seqSimStrInt.is_file(): 
        blast_interactome_protein_pairs(interactomeSequence, 
                                        structuralInteractome, 
                                        degreeStrInt, 
                                        IDMappingFile, 
                                        sga_GI_matrix, 
                                        blastTempFile, 
                                        querySeqTempFile, 
                                        subjectSeqTempFile, 
                                        seqSimStrInt)

    if not confounding_sequence_similarity_local_new.is_file(): 
        produce_sequence_similarity_covered_interfaces(structuralInteractomeWithInteractingResidues, 
                                                       degreeStrInt, 
                                                       IDMappingFile, 
                                                       sga_GI_matrix, 
                                                       protSeqDict, 
                                                       hubCutOff, 
                                                       TempFile_1, 
                                                       TempFile_2, 
                                                       blastTempFile, 
                                                       alignment_cutoff, 
                                                       confounding_sequence_similarity_local_new)

    # Structural similarity
    if not confounding_structural_similarity.is_file(): 
        produce_structural_similarity_whole_protein(structuralInteractome, 
                                                    degreeStrInt, 
                                                    IDMappingFile,
                                                    sga_GI_matrix,  
                                                    ppiModelDir, 
                                                    hubCutOff, 
                                                    TempFile_1, 
                                                    TempFile_2, 
                                                    TMalign_file, 
                                                    TMalignDir, 
                                                    confounding_structural_similarity)

    if not confounding_structural_similarity_local_new.is_file(): 
        produce_structural_similarity_extended_interfaces(structuralInteractomeWithInteractingResidues, 
                                                          degreeStrInt, 
                                                          IDMappingFile, 
                                                          sga_GI_matrix, 
                                                          protSeqDict, 
                                                          ppiModelDir, 
                                                          hubCutOff, 
                                                          TempFile_1, 
                                                          TempFile_2, 
                                                          TMalign_file, 
                                                          TMalignDir, 
                                                          strSim_cutoff, 
                                                          confounding_structural_similarity_local_new)

if __name__ == '__main__':
    main()