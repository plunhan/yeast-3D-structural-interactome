import os
import pandas as pd 
import numpy as np 
from pathlib import Path 
from analysis_tools import (count_degree_strInt,
                            generate_gene_list_strInt)

def main():

    # reference interactome name
    interactome_name = 'Yeast'

    # parent directory of all data files
    dataDir = Path('../data')

    # parent directory of all processed data files
    procDir = dataDir / 'processed'

    # parent directory of Aim 1 results
    rd1Dir = procDir / 'round_1'

    # parent directory of co-crystal analysis
    cocrystalDir = rd1Dir / 'cocrystal'

    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'validation'

    # input data files
    structuralInteractome = modelBasedDir / 'structural_interactome.txt'

    # output files
    degreeCount_strInt = cocrystalDir / 'degree_count_structural_interactome.txt' # to verify if PPI degree follows a power law distribution
    geneList = cocrystalDir / 'gene_list.txt' # to investigate if genes in the structural interactome are enriched for specific functions

    if not rd1Dir.exists():
        os.makedirs(rd1Dir)

    if not cocrystalDir.exists(): 
        os.makedirs(cocrystalDir)

    if not degreeCount_strInt.is_file():
        print('Counting PPI degree for genes in the structural interactome...')
        count_degree_strInt(structuralInteractome, degreeCount_strInt)

    if not geneList.is_file():
        print('Generating the list of genes within the structural interactome...')
        generate_gene_list_strInt(structuralInteractome, geneList)

if __name__ == '__main__':
    main()