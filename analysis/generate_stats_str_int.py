import os
import pandas as pd 
import numpy as np 
from pathlib import Path 
from analysis_tools import count_degree_strInt

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

    if not aim1Dir.exists():
        os.makedirs(aim1Dir)

    if not degreeCount_strInt.is_file():
        print('Counting PPI degree for genes in the structural interactome...')
        count_degree_strInt(structuralInteractome, degreeCount_strInt)

if __name__ == '__main__':
    main()