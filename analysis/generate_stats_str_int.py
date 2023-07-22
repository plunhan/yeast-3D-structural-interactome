import os
import pandas as pd 
import numpy as np 
from pathlib import Path 
from analysis_tools import (count_degree_refInt, 
                            count_degree_strInt, 
                            count_hub_for_each_number_of_interfaces, 
                            generate_gene_list_strInt, 
                            summarize_stats_strInt, 
                            degree_comparison)

def main():

    # reference interactome name
    interactome_name = 'Yeast'

    # parent directory of all data files
    dataDir = Path('../data')

    # parent directory of all processed data files
    procDir = dataDir / 'processed'

    # parent directory of Aim 1 results
    aim1Dir = procDir / 'aim_1'

    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'

    # input data files
    referenceInteractome = procDir / 'yeast_reference_interactome.txt'
    structuralInteractome = modelBasedDir / 'structural_interactome.txt'

    # output files
    degreeCount_strInt = aim1Dir / 'degree_count_structural_interactome.txt' # to verify if PPI degree follows a power law distribution
    degreeCount_refInt = aim1Dir / 'degree_count_reference_interactome.txt' # to verify if PPI degree follows a power law distribution
    degreeComparison = aim1Dir / 'degree_comparison.txt' # compare PPI degree between strInt and refInt for same protein
    geneList = aim1Dir / 'gene_list.txt' # to investigate if genes in the structural interactome are enriched for specific functions
    statsSummary = aim1Dir / 'stats_summary.txt' # to summarize stats derived from the structural interactome
    hubList = aim1Dir / 'hub_list.txt'

    if not aim1Dir.exists():
        os.makedirs(aim1Dir)

    if not degreeCount_refInt.is_file():
        print('Counting PPI degree for genes in the reference interactome...')
        count_degree_refInt(referenceInteractome, degreeCount_refInt)

    if not degreeCount_strInt.is_file():
        print('Counting PPI degree for genes in the structural interactome...')
        count_degree_strInt(structuralInteractome, degreeCount_strInt)

    if not degreeComparison.is_file():
        print('Generating file comparing PPI degree...')
        degree_comparison(degreeCount_refInt, degreeCount_strInt, degreeComparison)

    if not geneList.is_file():
        print('Generating the list of genes within the structural interactome...')
        generate_gene_list_strInt(structuralInteractome, geneList)

    if not statsSummary.is_file():
        print('Summarizing statistics for the structural interactome...')
        summarize_stats_strInt(structuralInteractome, referenceInteractome, statsSummary)

if __name__ == '__main__':
    main()