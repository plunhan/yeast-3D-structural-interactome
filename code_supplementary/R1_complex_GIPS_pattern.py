import os
from pathlib import Path
from import_tools import map_protein_id_to_locus_id
import pandas as pd
from scipy.stats import ranksums
import numpy as np

def calculate_GIPS(GI_profile, protein_1, protein_2): 
    if (protein_1 not in GI_profile.columns) or (protein_2 not in GI_profile.columns): 
        return False
    return GI_profile[[protein_1, protein_2]].corr().iloc[0,1]

def calculate_average_GIPS(protein, subunit_ls, GI_profile): 
    gips_ls = []
    for subunit in subunit_ls: 
        gips = calculate_GIPS(GI_profile, protein, subunit)
        if not gips or np.isnan(gips): 
            continue
        gips_ls.append(gips)
    if not len(gips_ls): 
        return -100
    return sum(gips_ls) / len(gips_ls)

def calculate_pairwise_average_GIPS(row, complex_dict, GI_profile): 
    hub = row['hub']
    partner_1 = row['partner_1']
    partner_2 = row['partner_2']
    subunit_ls = complex_dict[row['complex']]
    if len(subunit_ls) == 3: 
        return -100, -100
    subunit_ls_1 = [subunit for subunit in subunit_ls if subunit not in [hub, partner_2]]
    subunit_ls_2 = [subunit for subunit in subunit_ls if subunit not in [hub, partner_1]]
    average_gips_1 = calculate_average_GIPS(partner_1, subunit_ls_1, GI_profile)
    average_gips_2 = calculate_average_GIPS(partner_2, subunit_ls_2, GI_profile)
    if np.isnan(average_gips_1) or np.isnan(average_gips_2): 
        return -100, -100
    return average_gips_1, average_gips_2

def check_cocomplex(tuple_value, complex_dict):
    for complexName, protein_ls in complex_dict.items(): 
        if all(string in protein_ls for string in tuple_value): 
            return complexName
    return 'Non-cocomplex'

def create_tuple(row): 
    return (row['hub'], row['partner_1'], row['partner_2'])

def main(): 

    # Cut-off to define "large" interfacial overlap
    overlap_cutoff = 10

    # Cut-off to define "dissimilar" genetic interaction profiles
    GIPS_cutoff = 0.2

    interactome_name = 'Yeast'

    # parent directory of all data files
    dataDir = Path('../data')

    # parent directory of all external files
    extDir = dataDir / 'external'

    # parent directory of all processed data files
    procDir = dataDir / 'processed'

    # parent directory of round 1 results
    rd1Dir = procDir / 'round_1'

    # parent directory of Aim 2 results
    aim2Dir = procDir / 'aim_2'

    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based' 

    complexomeFile = extDir / '559292.tsv'
    complexome = pd.read_csv(complexomeFile, sep = '\t')
    complexome = complexome[['#Complex ac', 'Identifiers (and stoichiometry) of molecules in complex']]
    complexome.columns = ['Complex_name', 'Subunits']

    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)

    complex_dict = dict()
    list_of_lists = []

    for index, row in complexome.iterrows(): 
        complexName = row['Complex_name']
        proteins_stoc = row['Subunits'].split('|')
        proteins = [entry.split('(')[0] for entry in proteins_stoc]
        protein_ls = [IDMappingDict[protein] for protein in proteins if protein in IDMappingDict.keys()]
        protein_ls.sort()
        complex_dict[complexName] = protein_ls
        list_of_lists.append(protein_ls)

    InPathGI = procDir / 'sga_GI_rel_matrix.txt'
    GI_profile = pd.read_table(InPathGI)
    corResultFile = aim2Dir / 'overlapping_residue_count_vs_genetic_GIPS.txt'
    corResult = pd.read_table(corResultFile)
    corResult = corResult[corResult['ORC'] > overlap_cutoff][corResult['GIPS'] < GIPS_cutoff]
    corResult = corResult.dropna()
    corResult['index'] = corResult.apply(create_tuple, axis=1)
    corResult['complex'] = corResult['index'].apply(lambda x: check_cocomplex(x, complex_dict))
    corResult = corResult[corResult['complex'] != 'Non-cocomplex']
    corResult['complex_members'] = corResult['complex'].apply(lambda x: len(complex_dict[x]))
    corResult['GIPS_partner'] = corResult.apply(lambda row: calculate_pairwise_average_GIPS(row, complex_dict, GI_profile), axis=1)
    
    def split_tuple(row):
        return pd.Series(row)

    corResult[['GIPS_partner_1', 'GIPS_partner_2']] = corResult['GIPS_partner'].apply(split_tuple)
    corResult = corResult.drop(['index', 'GIPS_partner', 'PPI_degree'], axis=1)
    print(corResult)

    outPath = rd1Dir / 'R1_complex_GIPS_pattern.txt'
    corResult.to_csv(outPath, sep='\t')

if __name__ == '__main__':
    main()