import os
import pandas as pd 
import numpy as np 
from pathlib import Path 
from analysis_tools import count_degree_strInt, correlation_ORC_GIPS

def generate_chain_map_dict(strIntChainMap):
    chainMap_dict = dict(zip(zip(strIntChainMap['Query'], strIntChainMap['Subject']), strIntChainMap['Expect']))
    return chainMap_dict

def assign_quality(row, chainMap_dict, Evalue_cutoff): 
    alignment = row['Alignment_file_ID']
    p1, p2 = alignment.split('_')[0].split('=')
    template, ch1, ch2 = alignment.split('_')[1].split('-')
    if chainMap_dict[(p1, '_'.join([template, ch1.replace('!', '')]))] <= Evalue_cutoff and chainMap_dict[(p2, '_'.join([template, ch2.replace('!', '')]))] <= Evalue_cutoff: 
        return 'high'
    else: 
        return 'low'

def filter_high_quality_models(strInt, strIntChainMap, Evalue_cutoff):
    chainMap_dict = generate_chain_map_dict(strIntChainMap)
    strInt['quality'] = strInt.apply(lambda row: assign_quality(row, chainMap_dict, Evalue_cutoff), axis=1)
    return strInt[strInt['quality'] == 'high']

def main():

    Evalue_cutoff = 0

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

    # parent directory of Aim 1 results
    rd1Dir = procDir / 'round_1'

    # parent directory of co-crystal analysis
    cocrystalDir = rd1Dir / 'cocrystal'

    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name

    # directory of processed model-related data files specific to interactome
    validationModelBasedDir = interactomeDir / 'validation' / 'model_based'
    templateBasedDir = interactomeDir / 'template_based'
    modelBasedDir = interactomeDir / 'model_based'

    # input data files
    structuralInteractome_cocrystal = validationModelBasedDir / 'structural_interactome.txt'
    structuralInteractome = modelBasedDir / 'structural_interactome.txt'
    strucInteractomeChainMapFile = templateBasedDir / 'struc_interactome_chain_map.txt'
    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'
    geneticInteractionProfile_rel = procDir / 'sga_GI_rel_matrix.txt'

    # output data files
    highQualityStrInt = cocrystalDir / 'structural_interactome_high_quality.txt'
    degreeCountStrInt = cocrystalDir / 'degree_count_structural_interactome_high_quality.txt'
    corOverlapResCountGIPS = cocrystalDir / 'correlation_cocrystal_high_quality.txt'

    # correlation between overlapping residue count and genetic interaction profile similarity
    strIntChainMap = pd.read_table(strucInteractomeChainMapFile)
    strInt = pd.read_table(structuralInteractome)
    strInt = filter_high_quality_models(strInt, strIntChainMap, Evalue_cutoff)
    strInt.to_csv(highQualityStrInt, sep='\t', index=False)
    count_degree_strInt(highQualityStrInt, degreeCountStrInt)

    if not corOverlapResCountGIPS.is_file(): 
        correlation_ORC_GIPS(highQualityStrInt, 
                             degreeCountStrInt, 
                             IDMappingFile, 
                             geneticInteractionProfile_rel, 
                             hubCutOff, 
                             corOverlapResCountGIPS)

if __name__ == '__main__':
    main()