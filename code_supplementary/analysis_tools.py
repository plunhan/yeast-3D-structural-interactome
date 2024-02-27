import pandas as pd 
import numpy as np 
from scipy.stats import pearsonr, ranksums
from import_tools import (map_protein_id_to_locus_id, 
                          map_protein_id_to_gene_name, 
                          map_protein_id_to_ORF_name_pombe, 
                          import_structural_interactome, 
                          import_structural_interactome_with_interacting_residues, 
                          import_hub_proteins, 
                          import_human_hub_proteins, 
                          import_sequence_similarity, 
                          import_GI_profile, 
                          import_coexpression_data)
from interface_tools import (get_partner_list, 
                             calculate_degree_of_interface_overlap, 
                             calculate_GIPS, 
                             calculate_GIPS_human, 
                             calculate_GIPS_pombe)
from collections import Counter

######## Aim 1 ######## 
def count_degree_refInt(InPath, OutPath): 
    refInt = pd.read_table(InPath) 
    refInt['Protein_1'], refInt['Protein_2'] = refInt.min(axis=1), refInt.max(axis=1) 
    refInt = refInt.drop_duplicates() 
    d = dict(Counter(list(refInt['Protein_1'])+list(refInt['Protein_2']))) 
    df = pd.DataFrame(list(d.items()), columns=['Protein', 'Degree']) 
    df.to_csv(OutPath, sep='\t', index=False) 

def count_degree_strInt(InPath, OutPath): 
    strInt = pd.read_table(InPath, usecols=['Protein_1', 'Protein_2']) 
    strInt['Protein_1'], strInt['Protein_2'] = strInt.min(axis=1), strInt.max(axis=1) 
    strInt = strInt[strInt['Protein_1'] != strInt['Protein_2']] 
    strInt = strInt.drop_duplicates() 
    d = dict(Counter(list(strInt['Protein_1'])+list(strInt['Protein_2']))) 
    df = pd.DataFrame(list(d.items()), columns=['Protein', 'Degree']) 
    df.to_csv(OutPath, sep='\t', index=False) 

def degree_comparison(IntPathDegreeRef, IntPathDegreeStr, OutPath): 
    refIntDegree = pd.read_table(IntPathDegreeRef) 
    refIntDegree = refIntDegree.set_index('Protein') 
    strIntDegree = pd.read_table(IntPathDegreeStr) 
    strIntDegree.columns = ['Protein', 'Degree_strInt'] 
    strIntDegree = strIntDegree.set_index('Protein') 
    strIntDegree['Degree_refInt'] = refIntDegree['Degree'] 
    strIntDegree = strIntDegree.reset_index() 
    strIntDegree.to_csv(OutPath, sep='\t') 

def generate_gene_list_strInt(InPath, OutPath): 
    strInt = pd.read_table(InPath, usecols=['Protein_1', 'Protein_2']) 
    geneList = list(set(strInt['Protein_1']).union(set(strInt['Protein_2']))) 
    with open(OutPath, 'w') as f: 
        for gene in geneList: 
            f.write(gene + '\n') 

def summarize_stats_strInt(InPathStrInt, InPathRefInt, OutPath): 
    # 
    refInt = pd.read_table(InPathRefInt) 
    refInt['Protein_1'], refInt['Protein_2'] = refInt.min(axis=1), refInt.max(axis=1) 
    refInt = refInt.drop_duplicates() 
    ppiRefInt = len(refInt) 
    protRefInt = len(set(refInt['Protein_1']).union(set(refInt['Protein_2']))) 
    strInt = pd.read_table(InPathStrInt, usecols=['Protein_1', 'Protein_2']) 
    strInt = strInt.drop_duplicates() 
    ppiStrInt = len(strInt) 
    protStrInt = len(set(strInt['Protein_1']).union(set(strInt['Protein_2']))) 
    percentage = ppiStrInt / ppiRefInt 
    with open(OutPath, 'w') as f: 
        f.write('There are {} PPIs between {} proteins within the reference interactome.\n'.format(ppiRefInt, protRefInt)) 
        f.write('There are {} PPIs between {} proteins within the structural interactome.\n'.format(ppiStrInt, protStrInt)) 
        f.write('The structural interactome covers {:2.2%} PPIs within the reference interactome.\n'.format(percentage)) 

######## Aim 2 ######## 
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
    result_all_df.dropna(inplace=True)
    print()
    print(result_all_df)
    print(pearsonr(result_all_df['ORC'], result_all_df['GIPS']))
    result_all_df.to_csv(OutPath, sep='\t', index=False, na_rep='NA')

def correlation_ORC_GIPS_pombe(InPathStrInt, 
                               InPathDegreeStr, 
                               InPathIDMappingFile, 
                               InPathGIP, 
                               hubCutOff, 
                               OutPath): 
    strInt = pd.read_table(InPathStrInt)
    print(strInt)
    degreeStrInt = pd.read_table(InPathDegreeStr)
    hubProteins = list(degreeStrInt[degreeStrInt['Degree'] >= hubCutOff]['Protein'])
    print(hubProteins)
    GI_profile = pd.read_table(InPathGIP, index_col=0)
    print(GI_profile)
    IDMappingDict = map_protein_id_to_ORF_name_pombe(InPathIDMappingFile)
    result_all = []
    count = 0
    for hub in hubProteins: 
        print(hub)
        partnerList = get_partner_list(strInt, hub)
        print(partnerList)
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        for partner_1, partner_2 in [(a, b) for idx, a in enumerate(partnerList) for b in partnerList[idx + 1:]]: 
            index = '_'.join([min(partner_1, partner_2), max(partner_1, partner_2)])
            ORC, ORCR = calculate_degree_of_interface_overlap(strInt_temp, hub, partner_1, partner_2)
            print(ORC)
            GIPS = calculate_GIPS_pombe(GI_profile, IDMappingDict[partner_1], IDMappingDict[partner_2])
            print(GIPS)
            if (GIPS == False) or (GIPS == np.nan): 
                continue
            result_all.append((hub, partner_1, partner_2, ORC, GIPS))
        count += 1
        # print('\r{} of {} hub proteins are processed. '.format(count, len(hubProteins)), end='')
    result_all_df = pd.DataFrame(result_all, columns=['hub', 'partner_1', 'partner_2', 'ORC', 'GIPS'])
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
    
######## Aim 3 ########
def human_correlation_ORC_GIPS(InPathStrInt, 
                               InPathDegreeStr,
                               InPathIDMappingFile, 
                               humanPhenotypicProfile, 
                               hubCutOff, 
                               OutPath, 
                               Ghadie=False): 
    strInt = pd.read_table(InPathStrInt)
    if not Ghadie: 
        strInt = strInt[['Protein_1', 'Protein_2', 'Protein1_Interfacial_Residues_Uniprot', 'Protein2_Interfacial_Residues_Uniprot']]
        strInt['Interfaces'] = strInt[['Protein1_Interfacial_Residues_Uniprot', 'Protein2_Interfacial_Residues_Uniprot']].agg('+'.join, axis=1)
    IDMappingDict = map_protein_id_to_gene_name(InPathIDMappingFile)
    strInt['Protein_1'], strInt['Protein_2'] = strInt['Protein_1'].map(IDMappingDict), strInt['Protein_2'].map(IDMappingDict)
    hubProteins, hubDict = import_human_hub_proteins(InPathDegreeStr, InPathIDMappingFile, hubCutOff)
    GI_profile = pd.read_csv(humanPhenotypicProfile, index_col=0)
    GI_profile.columns = [s.split(' ')[0] for s in GI_profile.columns]
    result_all = []
    count = 0
    for hub in hubProteins: 
        print(hub)
        if hub in ['CRYAB', 'CRYAA', 'U2AF2']: 
            continue
        partnerList = get_partner_list(strInt, hub)
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        for partner_1, partner_2 in [(a, b) for idx, a in enumerate(partnerList) for b in partnerList[idx + 1:]]: 
            index = '_'.join([min(partner_1, partner_2), max(partner_1, partner_2)])
            ORC, ORCR = calculate_degree_of_interface_overlap(strInt_temp, hub, partner_1, partner_2)
            GIPS = calculate_GIPS_human(GI_profile, partner_1, partner_2)
            if (GIPS == False) or (GIPS == np.nan): 
                continue
            result_all.append((hub, partner_1, partner_2, ORC, GIPS, hubDict[hub]))
        count += 1
        print('\r{} of {} hub proteins are processed. '.format(count, len(hubProteins)), end='')
    result_all_df = pd.DataFrame(result_all, columns=['hub', 'partner_1', 'partner_2', 'ORC', 'GIPS', 'PPI_degree'])
    result_all_df.to_csv(OutPath, sep='\t', index=False, na_rep='NA')