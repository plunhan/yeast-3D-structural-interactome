'''
For yeast, more than 800 novel structures of complexes could be confidently predicted by AlphaFold2 (PMID: 34762488). 
The current work should leverage on this important dataset to further strengthen the coverage of its structural interactome.
'''

from pathlib import Path
import pandas as pd
import os
from analysis_tools import count_degree_strInt, correlation_ORC_GIPS

def merge_structural_interactome(strIntPath, alphaFoldPath, mergedStrIntPath): 
    strInt = pd.read_csv(strIntPath, sep='\t')
    alphaFold = pd.read_csv(alphaFoldPath, sep='\t', index_col=0)
    strInt['label'] = strInt.apply(lambda row: '_'.join(sorted([row['Protein_1'], row['Protein_2']])), axis=1)
    alphaFold['label'] = alphaFold.apply(lambda row: '_'.join(sorted([row['Protein_1'], row['Protein_2']])), axis=1)
    strInt = strInt[~strInt['label'].isin(alphaFold['label'])]
    mergedStrInt = pd.concat([strInt, alphaFold])
    mergedStrInt.to_csv(mergedStrIntPath, sep='\t')

def main(): 

    hubCutOff = 2

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

    # parent directory of round 1 review results
    round1Dir = procDir / 'round_1'

    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name

    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'

    # input files
    strInt = modelBasedDir / 'structural_interactome.txt'
    alphaFold = round1Dir / 'structural_interactome_alphafold.txt'
    sga_GI_matrix = procDir / 'sga_GI_rel_matrix.txt'
    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'

    # output files
    mergedStrInt = round1Dir / 'R1_structural_interactome_merged.txt'
    degreeCountStrInt = round1Dir / 'R1_degree_count.txt'
    alphaFoldResult = round1Dir / 'R1_rel_ORC_vs_GIPS.txt'

    if not round1Dir.is_dir(): 
        os.mkdir(round1Dir)

    if not mergedStrInt.is_file(): 
        merge_structural_interactome(strInt, alphaFold, mergedStrInt)

    if not degreeCountStrInt.is_file(): 
        count_degree_strInt(mergedStrInt, degreeCountStrInt)

    if not alphaFoldResult.is_file(): 
        correlation_ORC_GIPS(mergedStrInt, 
                             degreeCountStrInt, 
                             IDMappingFile, 
                             sga_GI_matrix, 
                             hubCutOff, 
                             alphaFoldResult)

if __name__ == '__main__':
    main()