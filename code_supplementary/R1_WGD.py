'''
Moreover, given the whole-genome duplication in S cerevisiae, 
it would be interesting to investigate how interacting pairs that are gene duplicates affect the result.
'''

from pathlib import Path
import pandas as pd
import os

def filter_WGD(allResult, WGDs, resultWGD, resultNonWGD, method='WGD'): 
    result_all = pd.read_csv(allResult, sep='\t')
    if method == 'WGD': 
        table_WGD = pd.read_excel(WGDs)
        table_WGD = table_WGD[['Gene 1', 'Gene 2']]
    elif method == 'All': 
        table_WGD = pd.read_excel(WGDs)
        table_WGD = table_WGD[['locus tag A', 'locus tag B']]
        table_WGD.columns = ['Gene 1', 'Gene 2']
    else: 
        print('Please choose a proper method. ')
    WGD_ls = set(list(table_WGD.apply(lambda row: '_'.join(sorted([row['Gene 1'], row['Gene 2']])), axis=1)))
    result_all['label'] = result_all.apply(lambda row: '_'.join(sorted([row['partner_1'], row['partner_2']])), axis=1)
    result_wgd = result_all[result_all['label'].isin(WGD_ls)]
    if not resultWGD.is_file(): 
        result_wgd.to_csv(resultWGD, sep='\t')
    result_nonwgd = result_all[~result_all['label'].isin(WGD_ls)]
    if not resultNonWGD.is_file(): 
        result_nonwgd.to_csv(resultNonWGD, sep='\t')

def main(): 

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
    degreeStrInt = aim1Dir / 'degree_count_structural_interactome.txt'
    allResult = aim2Dir / 'rel_overlapping_residue_count_vs_genetic_GIPS.txt'
    WGDs = extDir / 'Byrne and Wolfe 2005.xlsx'
    SSDs = extDir / 'Fares et al. 2013.xls'

    # output files
    resultWGD = round1Dir / 'result_WGD.txt'
    resultNonWGD = round1Dir / 'result_nonWGD.txt'
    resultSSD = round1Dir / 'result_SSD_WGD.txt'
    resultNonSSD = round1Dir / 'result_non_SSD_WGD.txt'

    if not round1Dir.is_dir(): 
        os.mkdir(round1Dir)

    if not resultWGD.is_file() or not resultNonWGD.is_file(): 
        filter_WGD(allResult, WGDs, resultWGD, resultNonWGD, method='WGD')

    if not resultSSD.is_file() or not resultNonSSD.is_file(): 
        filter_WGD(allResult, SSDs, resultSSD, resultNonSSD, method='All')

if __name__ == '__main__':
    main()