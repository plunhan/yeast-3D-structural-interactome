import pandas as pd 
import numpy as np 
from scipy.stats import pearsonr, ranksums
from import_tools import (map_protein_id_to_locus_id, 
                          import_structural_interactome, 
                          import_hub_proteins, 
                          import_GI_profile)
from interface_tools import (get_partner_list, 
                             calculate_degree_of_interface_overlap, 
                             calculate_GIPS)
from collections import Counter

def count_degree_strInt(InPath, OutPath): 
    '''

    This function generates the degree (the number of PPI a protein has) for each protein within the structural interactome. 

    Input: 

        InPath (Path): path to structural interactome file

    Output: 

        OutPath (Path): the degree for each protein within the structural interactome
    
    '''
    strInt = pd.read_table(InPath, usecols=['Protein_1', 'Protein_2']) 
    strInt['Protein_1'], strInt['Protein_2'] = strInt.min(axis=1), strInt.max(axis=1) 
    strInt = strInt[strInt['Protein_1'] != strInt['Protein_2']] 
    strInt = strInt.drop_duplicates() 
    d = dict(Counter(list(strInt['Protein_1'])+list(strInt['Protein_2']))) 
    df = pd.DataFrame(list(d.items()), columns=['Protein', 'Degree']) 
    df.to_csv(OutPath, sep='\t', index=False) 

def correlation_ORC_GIPS(InPathStrInt, 
                         InPathDegreeStr, 
                         InPathIDMappingFile, 
                         InPathGIP, 
                         hubCutOff, 
                         OutPath): 
    '''

    This function maps overlapping residue count (ORC) and genetic interaction profile similarity (GIPS) 
    to each interacting pairs. 

    Input: 
        
        InPathStrInt (Path): path to structural interactome file
        InPathDegreeStr (Path): path to PPI degree count of proteins within structural interactome
        IDMappingFile (Path): path to ID mapping file
        InPathGIP (Path): path to genetic interaction profile matrix
        hubCutOff (int): PPI degree above which a protein is defined as a hub

    Output: 

        DataFrame: interactor pairs with their ORC and GIPS
    
    '''

    strInt = import_structural_interactome(InPathStrInt, InPathIDMappingFile)
    hubProteins, hubDict = import_hub_proteins(InPathDegreeStr, InPathIDMappingFile, hubCutOff)
    GI_profile = import_GI_profile(InPathGIP)
    result_all = []
    count = 0
    for hub in hubProteins: 
        partnerList = get_partner_list(strInt, hub)
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        for partner_1, partner_2 in [(a, b) for idx, a in enumerate(partnerList) for b in partnerList[idx + 1:]]: 
            index = '_'.join([min(partner_1, partner_2), max(partner_1, partner_2)])
            ORC, ORCR = calculate_degree_of_interface_overlap(strInt_temp, hub, partner_1, partner_2)
            GIPS = calculate_GIPS(GI_profile, partner_1, partner_2)
            if (GIPS == False) or (GIPS == np.nan): 
                continue
            result_all.append((hub, partner_1, partner_2, ORC, ORCR, GIPS, hubDict[hub]))
        count += 1
        print('\r{} of {} hub proteins are processed. '.format(count, len(hubProteins)), end='')
    result_all_df = pd.DataFrame(result_all, columns=['hub', 'partner_1', 'partner_2', 'ORC', 'ORCR', 'GIPS', 'PPI_degree'])
    result_all_df.to_csv(OutPath, sep='\t', index=False, na_rep='NA')

def correlation_ORC_SemSim(InPathStrInt, 
                           InPathDegreeStr, 
                           InPathIDMappingFile, 
                           InPathGIP, 
                           InPathBP, 
                           InPathMF, 
                           InPathCC, 
                           hubCutOff, 
                           OutPath): 
    '''

    This function maps overlapping residue count (ORC) and genetic interaction profile similarity (GIPS) 
    to each interacting pairs. 

    Input: 
        
        InPathStrInt (Path): path to structural interactome file
        InPathDegreeStr (Path): path to PPI degree count of proteins within structural interactome
        IDMappingFile (Path): path to ID mapping file
        InPathGIP (Path): path to genetic interaction profile matrix
        hubCutOff (int): PPI degree above which a protein is defined as a hub

    Output: 

        DataFrame: interactor pairs with their ORC and GIPS
    
    '''

    strInt = import_structural_interactome(InPathStrInt, InPathIDMappingFile)
    hubProteins, hubDict = import_hub_proteins(InPathDegreeStr, InPathIDMappingFile, hubCutOff)
    IDMappingDict = map_protein_id_to_locus_id(InPathIDMappingFile)
    GI_profile = import_GI_profile(InPathGIP)

    BP = pd.read_table(InPathBP)
    BP['obj_1'], BP['obj_2'] = BP['obj_1'].map(IDMappingDict), BP['obj_2'].map(IDMappingDict)
    BP['obj_1'], BP['obj_2'] = BP[['obj_1', 'obj_2']].min(axis = 1), BP[['obj_1', 'obj_2']].max(axis = 1)
    BP['index'] = BP[['obj_1', 'obj_2']].agg('_'.join, axis = 1)
    BP_dict = dict(zip(BP['index'], BP['ss']))

    MF = pd.read_table(InPathMF)
    MF['obj_1'], MF['obj_2'] = MF['obj_1'].map(IDMappingDict), MF['obj_2'].map(IDMappingDict)
    MF['obj_1'], MF['obj_2'] = MF[['obj_1', 'obj_2']].min(axis = 1), MF[['obj_1', 'obj_2']].max(axis = 1)
    MF['index'] = MF[['obj_1', 'obj_2']].agg('_'.join, axis = 1)
    MF_dict = dict(zip(MF['index'], MF['ss']))

    CC = pd.read_table(InPathCC)
    CC['obj_1'], CC['obj_2'] = CC['obj_1'].map(IDMappingDict), CC['obj_2'].map(IDMappingDict)
    CC['obj_1'], CC['obj_2'] = CC[['obj_1', 'obj_2']].min(axis = 1), CC[['obj_1', 'obj_2']].max(axis = 1)
    CC['index'] = CC[['obj_1', 'obj_2']].agg('_'.join, axis = 1)
    CC_dict = dict(zip(CC['index'], CC['ss']))
    
    print('finish reading')
    result_all = []
    count = 0
    for hub in hubProteins: 
        partnerList = get_partner_list(strInt, hub)
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        for partner_1, partner_2 in [(a, b) for idx, a in enumerate(partnerList) for b in partnerList[idx + 1:]]: 
            index = '_'.join([min(partner_1, partner_2), max(partner_1, partner_2)])
            bp = BP_dict[index] if index in BP_dict.keys() else 'NA'
            mf = MF_dict[index] if index in MF_dict.keys() else 'NA'
            cc = CC_dict[index] if index in CC_dict.keys() else 'NA'
            ORC, ORCR = calculate_degree_of_interface_overlap(strInt_temp, hub, partner_1, partner_2)
            GIPS = calculate_GIPS(GI_profile, partner_1, partner_2)
            if (GIPS == False) or (GIPS == np.nan): 
                continue
            result_all.append((index, hub, partner_1, partner_2, ORC, ORCR, GIPS, bp, mf, cc, hubDict[hub]))
        count += 1
        print('\r{} of {} hub proteins are processed. '.format(count, len(hubProteins)), end='')
    result_all_df = pd.DataFrame(result_all, columns=['index', 'hub', 'partner_1', 'partner_2', 'ORC', 'ORCR', 'GIPS', 'BP', 'MF', 'cc', 'PPI_degree'])
    result_all_df.to_csv(OutPath, sep='\t', index=False, na_rep='NA')