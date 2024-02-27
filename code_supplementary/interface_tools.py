# For a pair of proteins that interact with the same target protein, this module finds the relationship between: 
# 1. The degree of interface overlap between the pair of proteins
# 2. The sequence/functional similarity for this pair of proteins, calculated by genetic interaction profile similarity/coexpression
import os
import pandas as pd 
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, binned_statistic, ranksums
import matplotlib.pyplot as plt
from import_tools import map_protein_id_to_locus_id

def calculate_degree_of_interface_overlap(strInt, hub, protein_1, protein_2): 
    # Measuring "degree of interface overlap": 
    # 1. The number of overlapping residues; 
    # 2. The ratio of the number of interfacial residues, to the total number of interfacial residues targeted by either of the two interacting partners;
    # 3. The distance between the geometric center of the two interfaces. 
    strInt_1 = strInt[(strInt['Protein_1'] == protein_1) | (strInt['Protein_2'] == protein_1)]
    if strInt_1['Protein_1'].iloc[0] == hub: 
        interface_1 = set(strInt_1['Interfaces'].iloc[0].split('+')[0].split(','))
    else: 
        interface_1 = set(strInt_1['Interfaces'].iloc[0].split('+')[1].split(','))
    strInt_2 = strInt[(strInt['Protein_1'] == protein_2) | (strInt['Protein_2'] == protein_2)]
    if strInt_2['Protein_1'].iloc[0] == hub: 
        interface_2 = set(strInt_2['Interfaces'].iloc[0].split('+')[0].split(','))
    else: 
        interface_2 = set(strInt_2['Interfaces'].iloc[0].split('+')[1].split(','))
    return len(interface_1.intersection(interface_2)), len(interface_1.intersection(interface_2))/(len(interface_1.union(interface_2)))

def calculate_GIPS(GI_profile, protein_1, protein_2): 
    if (protein_1 not in GI_profile.columns) or (protein_2 not in GI_profile.columns): 
        return False
    return GI_profile[[protein_1, protein_2]].corr().iloc[0,1]

def calculate_GIPS_human(GI_profile, protein_1, protein_2): 
    if (protein_1 not in GI_profile.columns) or (protein_2 not in GI_profile.columns): 
        return False
    return pearsonr(GI_profile[protein_1], GI_profile[protein_2])[0]

def calculate_GIPS_pombe(GI_profile, protein_1, protein_2): 
    def in_column(GI_profile, protein): 
        return protein in GI_profile.columns
    def in_row(GI_profile, protein): 
        return protein in GI_profile.index
    if in_row(GI_profile, protein_1) and in_row(GI_profile, protein_2): 
        return GI_profile.loc[[protein_1, protein_2]].T.corr().iloc[0, 1]
    elif in_column(GI_profile, protein_1) and in_column(GI_profile, protein_2): 
        return GI_profile[[protein_1, protein_2]].corr().iloc[0, 1]
    return False

def calculate_GI_share(SGA_matrix, protein_1, protein_2): 
    if (protein_1 not in SGA_matrix.columns) or (protein_2 not in SGA_matrix.columns): 
        return False
    temp =  SGA_matrix[protein_1] * SGA_matrix[protein_2]
    return temp.sum()

def import_old_coexpression_data(InPath): 
    # Import co-expression data
    cx = pd.read_table(InPath)
    cx = cx.drop(0)
    cx = cx.drop(['NAME', 'GWEIGHT'], axis=1)
    cx = cx.drop_duplicates(subset=['YORF'], keep='first')
    cx = cx.set_index('YORF')
    return cx

def calculate_coexpression(coexpressionFile, protein_1, protein_2): 
    if (protein_1 not in coexpressionFile.index) or (protein_2 not in coexpressionFile.index): 
        return False
    coexpressionScore = pearsonr(coexpressionFile.loc[protein_1], coexpressionFile.loc[protein_2])[0]
    return coexpressionScore

def calculate_semantic_similarity(semsimDict, protein_1, protein_2): 
    protein_1, protein_2 = min([protein_1, protein_2]), max([protein_1, protein_2])
    index = '_'.join([protein_1, protein_2])
    if index not in semsimDict.keys(): 
        return False
    return semsimDict[index]

def calculate_sequence_similarity(seqsimDict, protein_1, protein_2): 
    protein_1, protein_2 = min([protein_1, protein_2]), max([protein_1, protein_2])
    index = '_'.join([protein_1, protein_2])
    if index not in seqsimDict.keys(): 
        return False
    return seqsimDict[index]

def get_partner_list(strInt, hub): 
    strInt_hub = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
    partnerList = list(set(strInt_hub['Protein_1']).union(set(strInt_hub['Protein_2'])))
    partnerList.remove(hub)
    partnerList.sort()
    return partnerList

def print_statistical_difference(result_nonparalogs_df, 
                                 result_paralogs_df, 
                                 result_WGD_df, 
                                 result_SSD_df): 
    print('degree_1')
    print(result_nonparalogs_df['degree_1'].mean(), result_paralogs_df['degree_1'].mean())
    print(ranksums(result_nonparalogs_df['degree_1'], result_paralogs_df['degree_1']))
    print(result_WGD_df['degree_1'].mean(), result_nonparalogs_df['degree_1'].mean())
    print(ranksums(result_WGD_df['degree_1'], result_nonparalogs_df['degree_1']))
    print(result_SSD_df['degree_1'].mean(), result_nonparalogs_df['degree_1'].mean())
    print(ranksums(result_SSD_df['degree_1'], result_nonparalogs_df['degree_1']))
    print(result_WGD_df['degree_1'].mean(), result_SSD_df['degree_1'].mean())
    print(ranksums(result_WGD_df['degree_1'], result_SSD_df['degree_1']))
    print('degree_2')
    print(result_nonparalogs_df['degree_2'].mean(), result_paralogs_df['degree_2'].mean())
    print(ranksums(result_nonparalogs_df['degree_2'], result_paralogs_df['degree_2']))
    print(result_WGD_df['degree_2'].mean(), result_nonparalogs_df['degree_2'].mean())
    print(ranksums(result_WGD_df['degree_2'], result_nonparalogs_df['degree_2']))
    print(result_SSD_df['degree_2'].mean(), result_nonparalogs_df['degree_2'].mean())
    print(ranksums(result_SSD_df['degree_2'], result_nonparalogs_df['degree_2']))
    print(result_WGD_df['degree_2'].mean(), result_SSD_df['degree_2'].mean())
    print(ranksums(result_WGD_df['degree_2'], result_SSD_df['degree_2']))
    print('coexpression')
    print(result_nonparalogs_df['coexpression'].mean(), result_paralogs_df['coexpression'].mean())
    print(ranksums(result_nonparalogs_df['coexpression'], result_paralogs_df['coexpression']))
    print(result_WGD_df['coexpression'].mean(), result_nonparalogs_df['coexpression'].mean())
    print(ranksums(result_WGD_df['coexpression'], result_nonparalogs_df['coexpression']))
    print(result_SSD_df['coexpression'].mean(), result_nonparalogs_df['coexpression'].mean())
    print(ranksums(result_SSD_df['coexpression'], result_nonparalogs_df['coexpression']))
    print(result_WGD_df['coexpression'].mean(), result_SSD_df['coexpression'].mean())
    print(ranksums(result_WGD_df['coexpression'], result_SSD_df['coexpression']))
    print('GIPS')
    print(result_nonparalogs_df['GIPS'].mean(), result_paralogs_df['GIPS'].mean())
    print(ranksums(result_nonparalogs_df['GIPS'], result_paralogs_df['GIPS']))
    print(result_WGD_df['GIPS'].mean(), result_nonparalogs_df['GIPS'].mean())
    print(ranksums(result_WGD_df['GIPS'], result_nonparalogs_df['GIPS']))
    print(result_SSD_df['GIPS'].mean(), result_nonparalogs_df['GIPS'].mean())
    print(ranksums(result_SSD_df['GIPS'], result_nonparalogs_df['GIPS']))
    print(result_WGD_df['GIPS'].mean(), result_SSD_df['GIPS'].mean())
    print(ranksums(result_WGD_df['GIPS'], result_SSD_df['GIPS']))
    print('SharedGIs')
    print(result_nonparalogs_df['SharedGIs'].mean(), result_paralogs_df['SharedGIs'].mean())
    print(ranksums(result_nonparalogs_df['SharedGIs'], result_paralogs_df['SharedGIs']))
    print(result_WGD_df['SharedGIs'].mean(), result_nonparalogs_df['SharedGIs'].mean())
    print(ranksums(result_WGD_df['SharedGIs'], result_nonparalogs_df['SharedGIs']))
    print(result_SSD_df['SharedGIs'].mean(), result_nonparalogs_df['SharedGIs'].mean())
    print(ranksums(result_SSD_df['SharedGIs'], result_nonparalogs_df['SharedGIs']))
    print(result_WGD_df['SharedGIs'].mean(), result_SSD_df['SharedGIs'].mean())
    print(ranksums(result_WGD_df['SharedGIs'], result_SSD_df['SharedGIs']))
    print('SemSim_mf')
    print(result_nonparalogs_df['SemSim_mf'].mean(), result_paralogs_df['SemSim_mf'].mean())
    print(ranksums(result_nonparalogs_df['SemSim_mf'], result_paralogs_df['SemSim_mf']))
    print(result_WGD_df['SemSim_mf'].mean(), result_nonparalogs_df['SemSim_mf'].mean())
    print(ranksums(result_WGD_df['SemSim_mf'], result_nonparalogs_df['SemSim_mf']))
    print(result_SSD_df['SemSim_mf'].mean(), result_nonparalogs_df['SemSim_mf'].mean())
    print(ranksums(result_SSD_df['SemSim_mf'], result_nonparalogs_df['SemSim_mf']))
    print(result_WGD_df['SemSim_mf'].mean(), result_SSD_df['SemSim_mf'].mean())
    print(ranksums(result_WGD_df['SemSim_mf'], result_SSD_df['SemSim_mf']))
    print('SemSim_bp')
    print(result_nonparalogs_df['SemSim_bp'].mean(), result_paralogs_df['SemSim_bp'].mean())
    print(ranksums(result_nonparalogs_df['SemSim_bp'], result_paralogs_df['SemSim_bp']))
    print(result_WGD_df['SemSim_bp'].mean(), result_nonparalogs_df['SemSim_bp'].mean())
    print(ranksums(result_WGD_df['SemSim_bp'], result_nonparalogs_df['SemSim_bp']))
    print(result_SSD_df['SemSim_bp'].mean(), result_nonparalogs_df['SemSim_bp'].mean())
    print(ranksums(result_SSD_df['SemSim_bp'], result_nonparalogs_df['SemSim_bp']))
    print(result_WGD_df['SemSim_bp'].mean(), result_SSD_df['SemSim_bp'].mean())
    print(ranksums(result_WGD_df['SemSim_bp'], result_SSD_df['SemSim_bp']))
    print('Identities')
    print(result_nonparalogs_df['Identities'].mean(), result_paralogs_df['Identities'].mean())
    print(ranksums(result_nonparalogs_df['Identities'], result_paralogs_df['Identities']))
    print(result_WGD_df['Identities'].mean(), result_nonparalogs_df['Identities'].mean())
    print(ranksums(result_WGD_df['Identities'], result_nonparalogs_df['Identities']))
    print(result_SSD_df['Identities'].mean(), result_nonparalogs_df['Identities'].mean())
    print(ranksums(result_SSD_df['Identities'], result_nonparalogs_df['Identities']))
    print(result_WGD_df['Identities'].mean(), result_SSD_df['Identities'].mean())
    print(ranksums(result_WGD_df['Identities'], result_SSD_df['Identities']))
    print('Positives')
    print(result_nonparalogs_df['Positives'].mean(), result_paralogs_df['Positives'].mean())
    print(ranksums(result_nonparalogs_df['Positives'], result_paralogs_df['Positives']))
    print(result_WGD_df['Positives'].mean(), result_nonparalogs_df['Positives'].mean())
    print(ranksums(result_WGD_df['Positives'], result_nonparalogs_df['Positives']))
    print(result_SSD_df['Positives'].mean(), result_nonparalogs_df['Positives'].mean())
    print(ranksums(result_SSD_df['Positives'], result_nonparalogs_df['Positives']))
    print(result_WGD_df['Positives'].mean(), result_SSD_df['Positives'].mean())
    print(ranksums(result_WGD_df['Positives'], result_SSD_df['Positives']))

def degree_of_interface_overlap_vs_functional_similarity(InPathStrInt, 
                                                         InPathDegree, 
                                                         InPathCoexpression, 
                                                         InPathGI, 
                                                         InPathSGAMatrix, 
                                                         InPathSemSim_mf, 
                                                         InPathSemSim_bp, 
                                                         InPathSeqSim, 
                                                         InPathAllParalogs, 
                                                         InPathWGD, 
                                                         IDMappingFile, 
                                                         OutPath, 
                                                         OutPathParalogs, 
                                                         OutPathWGD, 
                                                         OutPathSSD, 
                                                         OutPathNonparalogs): 
    # Import files
    strInt = pd.read_table(InPathStrInt)
    degreeStrInt = pd.read_table(InPathDegree)
    hubProteins = list(set(degreeStrInt[degreeStrInt['Degree'] >= 4]['Protein']))
    coexpressionFile = import_old_coexpression_data(InPathCoexpression)
    GI_profile = pd.read_table(InPathGI)
    SGA_matrix = pd.read_table(InPathSGAMatrix, index_col=0)
    semsim_mf = pd.read_table(InPathSemSim_mf)
    semsim_bp = pd.read_table(InPathSemSim_bp)
    seqsim = pd.read_table(InPathSeqSim)
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    hubProteins = [IDMappingDict[x] for x in hubProteins if x in IDMappingDict.keys()]

    # Process structural interactome
    strInt['Protein_1'] = strInt['Protein_1'].map(IDMappingDict)
    strInt['Protein_2'] = strInt['Protein_2'].map(IDMappingDict)
    strInt = strInt[strInt['Protein_1'] != strInt['Protein_2']]
    strInt = strInt.drop_duplicates()

    # Process semantic similarity file
    semsim_mf['obj_1'] = semsim_mf['obj_1'].map(IDMappingDict)
    semsim_mf['obj_2'] = semsim_mf['obj_2'].map(IDMappingDict)
    semsim_mf['obj_1'], semsim_mf['obj_2'] = semsim_mf[['obj_1', 'obj_2']].min(axis=1), semsim_mf[['obj_1', 'obj_2']].max(axis=1)
    semsim_mf['index'] = semsim_mf[['obj_1', 'obj_2']].agg('_'.join, axis=1)
    semsim_mf = semsim_mf.drop_duplicates()
    semsim_mf_Dict = dict([(key, value) for key, value in zip(semsim_mf['index'], semsim_mf['ss'])])
    semsim_bp['obj_1'] = semsim_bp['obj_1'].map(IDMappingDict)
    semsim_bp['obj_2'] = semsim_bp['obj_2'].map(IDMappingDict)
    semsim_bp['obj_1'], semsim_bp['obj_2'] = semsim_bp[['obj_1', 'obj_2']].min(axis=1), semsim_bp[['obj_1', 'obj_2']].max(axis=1)
    semsim_bp['index'] = semsim_bp[['obj_1', 'obj_2']].agg('_'.join, axis=1)
    semsim_bp = semsim_bp.drop_duplicates()
    semsim_bp_Dict = dict([(key, value) for key, value in zip(semsim_bp['index'], semsim_bp['ss'])])

    # Process sequence similarity file
    seqsim['Protein_1'] = seqsim['Protein_1'].map(IDMappingDict)
    seqsim['Protein_2'] = seqsim['Protein_2'].map(IDMappingDict)
    seqsim['Protein_1'], seqsim['Protein_2'] = seqsim[['Protein_1', 'Protein_2']].min(axis=1), seqsim[['Protein_1', 'Protein_2']].max(axis=1)
    seqsim['index'] = seqsim[['Protein_1', 'Protein_2']].agg('_'.join, axis=1)
    identitiesDict = dict([(key, value) for key, value in zip(seqsim['index'], seqsim['Identities'])])
    positivesDict = dict([(key, value) for key, value in zip(seqsim['index'], seqsim['Positives'])])

    # Generate paralogs within structural interactome
    WGD_list, SSD_list = check_paralogs(strInt, degreeStrInt, hubProteins, IDMappingDict, InPathAllParalogs, InPathWGD)

    result_all = []
    result_paralogs = []
    result_nonparalogs = []
    result_WGD = []
    result_SSD = []
    count = 0
    for hub in hubProteins: 
        partnerList = get_partner_list(strInt, hub)
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        for i in range(len(partnerList)-1): 
            for j in range(i, len(partnerList)): 
                if partnerList[i] == partnerList[j]: 
                    continue
                index = '_'.join([min(partnerList[i], partnerList[j]), max(partnerList[i], partnerList[j])])
                degree_1, degree_2 = calculate_degree_of_interface_overlap(strInt_temp, hub, partnerList[i], partnerList[j])
                coexpressionScore = calculate_coexpression(coexpressionFile, partnerList[i], partnerList[j])
                if coexpressionScore == False: 
                    continue
                GIPS = calculate_GIPS(GI_profile, partnerList[i], partnerList[j])
                if GIPS == False: 
                    continue
                GI_share = calculate_GI_share(SGA_matrix, partnerList[i], partnerList[j])
                if GI_share == False: 
                    continue
                SS_mf = calculate_semantic_similarity(semsim_mf_Dict, partnerList[i], partnerList[j])
                if SS_mf == False: 
                    continue
                SS_bp = calculate_semantic_similarity(semsim_bp_Dict, partnerList[i], partnerList[j])
                if SS_bp == False: 
                    continue
                identities = calculate_sequence_similarity(identitiesDict, partnerList[i], partnerList[j])
                if identities == False: 
                    continue
                positives = calculate_sequence_similarity(positivesDict, partnerList[i], partnerList[j])
                if positives == False: 
                    continue
                result_all.append((hub, partnerList[i], partnerList[j], degree_1, degree_2, coexpressionScore, GIPS, GI_share, SS_mf, SS_bp, identities, positives))
                if index in set(WGD_list): 
                    result_paralogs.append((hub, partnerList[i], partnerList[j], degree_1, degree_2, coexpressionScore, GIPS, GI_share, SS_mf, SS_bp, identities, positives, 'WGD'))
                    result_WGD.append((hub, partnerList[i], partnerList[j], degree_1, degree_2, coexpressionScore, GIPS, GI_share, SS_mf, SS_bp, identities, positives))
                elif index in set(SSD_list): 
                    result_paralogs.append((hub, partnerList[i], partnerList[j], degree_1, degree_2, coexpressionScore, GIPS, GI_share, SS_mf, SS_bp, identities, positives, 'SSD'))
                    result_SSD.append((hub, partnerList[i], partnerList[j], degree_1, degree_2, coexpressionScore, GIPS, GI_share, SS_mf, SS_bp, identities, positives))
                else: 
                    result_nonparalogs.append((hub, partnerList[i], partnerList[j], degree_1, degree_2, coexpressionScore, GIPS, GI_share, SS_mf, SS_bp, identities, positives))
        count += 1
        # print("{} of {} hub proteins have been processed...".format(count, len(hubProteins)))
    result_df = pd.DataFrame(result_all, columns=['hub', 'partner_1', 'partner_2', 'degree_1', 'degree_2', 'coexpression', 'GIPS', 'SharedGIs', 'SemSim_mf', 'SemSim_bp', 'Identities', 'Positives'])
    result_paralogs_df = pd.DataFrame(result_paralogs, columns=['hub', 'partner_1', 'partner_2', 'degree_1', 'degree_2', 'coexpression', 'GIPS', 'SharedGIs', 'SemSim_mf', 'SemSim_bp', 'Identities', 'Positives', 'Paralog type'])
    result_WGD_df = pd.DataFrame(result_WGD, columns=['hub', 'partner_1', 'partner_2', 'degree_1', 'degree_2', 'coexpression', 'GIPS', 'SharedGIs', 'SemSim_mf', 'SemSim_bp', 'Identities', 'Positives'])
    result_SSD_df = pd.DataFrame(result_SSD, columns=['hub', 'partner_1', 'partner_2', 'degree_1', 'degree_2', 'coexpression', 'GIPS', 'SharedGIs', 'SemSim_mf', 'SemSim_bp', 'Identities', 'Positives'])
    result_nonparalogs_df = pd.DataFrame(result_nonparalogs, columns=['hub', 'partner_1', 'partner_2', 'degree_1', 'degree_2', 'coexpression', 'GIPS', 'SharedGIs', 'SemSim_mf', 'SemSim_bp', 'Identities', 'Positives'])
    print_statistical_difference(result_nonparalogs_df, result_paralogs_df, result_WGD_df, result_SSD_df)
    result_WGD_df.to_csv(OutPathWGD, sep='\t', index=False)
    result_SSD_df.to_csv(OutPathSSD, sep='\t', index=False)
    result_nonparalogs_df.to_csv(OutPathNonparalogs, sep='\t', index=False)
    if not OutPath.is_file(): 
        result_df.to_csv(OutPath, sep='\t', index=False)
    if not OutPathParalogs.is_file():
        result_paralogs_df.to_csv(OutPathParalogs, sep='\t', index=False)

def binning_data(InPath, OutPath): 
    result_df = pd.read_table(InPath)
    result_df = result_df.dropna()
    upperright = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 27, 29, 31, 33, 34, 36, 39, 45, 67]
    middle = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 28, 30, 32, 33, 35, 37.5, 42, 56]
    bin_means_coexpression, bin_edges, binnumber = binned_statistic(result_df['degree_1'], result_df['coexpression'], 'mean', bins = upperright)
    bin_means_GIPS, bin_edges, binnumber = binned_statistic(result_df['degree_1'], result_df['GIPS'], 'mean', bins = upperright)
    bin_means_SemSim, bin_edges, binnumber = binned_statistic(result_df['degree_1'], result_df['SemSim'], 'mean', bins = upperright)
    bin_means_Positives, bin_edges, binnumber = binned_statistic(result_df['degree_1'], result_df['Positives'], 'mean', bins = upperright)
    bin_means_Identities, bin_edges, binnumber = binned_statistic(result_df['degree_1'], result_df['Identities'], 'mean', bins = upperright)
    output = pd.DataFrame({'edges': middle, 
                           'coexpression': bin_means_coexpression, 
                           'GIPS': bin_means_GIPS, 
                           'SemSim': bin_means_SemSim, 
                           'Identities': bin_means_Identities, 
                           'Positives': bin_means_Positives}, 
                           columns = ['edges', 'coexpression', 'GIPS', 'SemSim', 'Identities', 'Positives'])
    output.to_csv(OutPath, sep = '\t')

def get_all_partner_list(strInt, hubProteins): 
    allPartnerList = []
    for hub in hubProteins: 
        partnerList = get_partner_list(strInt, hub)
        if len(partnerList) <= 1: 
            continue
        for i in range(len(partnerList)-1): 
            for j in range(i+1, len(partnerList)): 
                if partnerList[i] == partnerList[j]: 
                    continue
                allPartnerList.append('_'.join([min(partnerList[i], partnerList[j]), max(partnerList[i], partnerList[j])]))
    return allPartnerList


def check_paralogs(strInt, degreeStrInt, hubProteins, IDMappingDict, InPathAllParalogs, InPathWGD): 
    # import files
    allParalogs = pd.read_excel(InPathAllParalogs)
    allParalogs = allParalogs[['locus tag A', 'locus tag B']]
    WGD = pd.read_excel(InPathWGD)
    WGD = WGD[['Gene 1', 'Gene 2']]
    allPartnerList = get_all_partner_list(strInt, hubProteins)
    print('There are {} pairs of interactors of hub proteins in the structural interactome. '.format(len(set(allPartnerList))))
    allParalogs['locus tag A'], allParalogs['locus tag B'] = allParalogs[['locus tag A', 'locus tag B']].min(axis=1), allParalogs[['locus tag A', 'locus tag B']].max(axis=1)
    allParalogs['index'] = allParalogs.agg('_'.join, axis=1)
    WGD['Gene 1'], WGD['Gene 2'] = WGD[['Gene 1', 'Gene 2']].min(axis=1), WGD[['Gene 1', 'Gene 2']].max(axis=1)
    WGD['index'] = WGD.agg('_'.join, axis=1)
    WGD_list = list(set(WGD['index']))
    SSD_list = list(set(allParalogs['index']).difference(set(WGD_list)))
    print('There are {} WGD paralog pairs and {} SSD paralog pairs. '.format(len(WGD_list), len(SSD_list)))
    print('There are {} WGD paralog pairs and {} SSD paralog pairs in the structural interactome. '.format(len(set(WGD_list).intersection(set(allPartnerList))), len(set(SSD_list).intersection(set(allPartnerList)))))
    return WGD_list, SSD_list

def generate_partner_evalue_distribution(InPathStrInt, 
                                         InPathDegree, 
                                         InPathAllParalogs, 
                                         IDMappingFile, 
                                         InPathEvalue, 
                                         OutPathPartners, 
                                         OutPathParalogs,  
                                         OutPathNonparalogs): 
    # Import files
    strInt = pd.read_table(InPathStrInt)
    degreeStrInt = pd.read_table(InPathDegree)
    hubProteins = list(set(degreeStrInt[degreeStrInt['Degree'] >= 4]['Protein']))
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    hubProteins = [IDMappingDict[x] for x in hubProteins if x in IDMappingDict.keys()]
    Evalue = pd.read_table(InPathEvalue)
    allParalogs = pd.read_excel(InPathAllParalogs)
    allParalogs = allParalogs[['locus tag A', 'locus tag B']]

    # Process structural interactome
    strInt['Protein_1'] = strInt['Protein_1'].map(IDMappingDict)
    strInt['Protein_2'] = strInt['Protein_2'].map(IDMappingDict)
    strInt = strInt[strInt['Protein_1'] != strInt['Protein_2']]
    strInt = strInt.drop_duplicates()

    # Process e-value file
    Evalue['Protein_1'] = Evalue['Protein_1'].map(IDMappingDict)
    Evalue['Protein_2'] = Evalue['Protein_2'].map(IDMappingDict)
    Evalue['Protein_1'], Evalue['Protein_2'] = Evalue[['Protein_1', 'Protein_2']].min(axis=1), Evalue[['Protein_1', 'Protein_2']].max(axis=1)
    Evalue['index'] = Evalue[['Protein_1', 'Protein_2']].agg('_'.join, axis=1)

    # Process paralog file
    allParalogs['locus tag A'], allParalogs['locus tag B'] = allParalogs[['locus tag A', 'locus tag B']].min(axis=1), allParalogs[['locus tag A', 'locus tag B']].max(axis=1)
    allParalogs['index'] = allParalogs.agg('_'.join, axis=1)

    # Get all interactor pairs
    allPartnerSet = set()
    paralogPartnerSet = set(allParalogs['index'])
    for hub in hubProteins: 
        partnerList = get_partner_list(strInt, hub)
        for i in range(len(partnerList)-1): 
            for j in range(i+1, len(partnerList)): 
                if partnerList[i] == partnerList[j]: 
                    continue
                else: 
                    allPartnerSet.add('_'.join([min(partnerList[i], partnerList[j]), max(partnerList[i], partnerList[j])]))
    print(allPartnerSet)
    print(Evalue['index'])
    Evalue_partner = Evalue[Evalue['index'].isin(allPartnerSet)]
    Evalue_partner.to_csv(OutPathPartners, sep='\t')
    Evalue_paralog = Evalue[Evalue['index'].isin(paralogPartnerSet)]
    Evalue_paralog.to_csv(OutPathParalogs, sep='\t')
    Evalue_nonparalog = Evalue[~Evalue['index'].isin(paralogPartnerSet)]
    Evalue_nonparalog.to_csv(OutPathNonparalogs, sep='\t')

def calculate_evalue(DataFrame, protein_1, protein_2): 
    index = '_'.join([min(protein_1, protein_2), max(protein_1, protein_2)])
    if index in set(DataFrame['index']): 
        return DataFrame[DataFrame['index'] == index]['Evalue'].iloc[0]
    else: 
        return False

def degree_of_interface_overlap_vs_functional_similarity_new(InPathStrInt, 
                                                             InPathDegree, 
                                                             InPathCoexpression, 
                                                             InPathGI, 
                                                             InPathSGAMatrix, 
                                                             InPathSemSim_mf, 
                                                             InPathSemSim_bp, 
                                                             InPathSeqSim, 
                                                             IDMappingFile, 
                                                             InPathDistEvalueAll, 
                                                             InPathDistEvalueParalogs, 
                                                             InPathDistEvalueNonparalogs, 
                                                             OutPathAll, 
                                                             OutPathParalogs, 
                                                             OutPathNonparalogs):
    # Import files
    strInt = pd.read_table(InPathStrInt)
    degreeStrInt = pd.read_table(InPathDegree)
    hubProteins = list(set(degreeStrInt[degreeStrInt['Degree'] >= 4]['Protein']))
    coexpressionFile = import_old_coexpression_data(InPathCoexpression)
    GI_profile = pd.read_table(InPathGI)
    SGA_matrix = pd.read_table(InPathSGAMatrix, index_col=0)
    semsim_mf = pd.read_table(InPathSemSim_mf)
    semsim_bp = pd.read_table(InPathSemSim_bp)
    seqsim = pd.read_table(InPathSeqSim)
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    hubProteins = [IDMappingDict[x] for x in hubProteins if x in IDMappingDict.keys()]
    distEvalue_all = pd.read_table(InPathDistEvalueAll)
    distEvalue_paralogs = pd.read_table(InPathDistEvalueParalogs)
    distEvalue_nonparalogs = pd.read_table(InPathDistEvalueNonparalogs)

    # Process structural interactome
    strInt['Protein_1'] = strInt['Protein_1'].map(IDMappingDict)
    strInt['Protein_2'] = strInt['Protein_2'].map(IDMappingDict)
    strInt = strInt[strInt['Protein_1'] != strInt['Protein_2']]
    strInt = strInt.drop_duplicates()

    # Process semantic similarity file
    semsim_mf['obj_1'] = semsim_mf['obj_1'].map(IDMappingDict)
    semsim_mf['obj_2'] = semsim_mf['obj_2'].map(IDMappingDict)
    semsim_mf['obj_1'], semsim_mf['obj_2'] = semsim_mf[['obj_1', 'obj_2']].min(axis=1), semsim_mf[['obj_1', 'obj_2']].max(axis=1)
    semsim_mf['index'] = semsim_mf[['obj_1', 'obj_2']].agg('_'.join, axis=1)
    semsim_mf = semsim_mf.drop_duplicates()
    semsim_mf_Dict = dict([(key, value) for key, value in zip(semsim_mf['index'], semsim_mf['ss'])])
    semsim_bp['obj_1'] = semsim_bp['obj_1'].map(IDMappingDict)
    semsim_bp['obj_2'] = semsim_bp['obj_2'].map(IDMappingDict)
    semsim_bp['obj_1'], semsim_bp['obj_2'] = semsim_bp[['obj_1', 'obj_2']].min(axis=1), semsim_bp[['obj_1', 'obj_2']].max(axis=1)
    semsim_bp['index'] = semsim_bp[['obj_1', 'obj_2']].agg('_'.join, axis=1)
    semsim_bp = semsim_bp.drop_duplicates()
    semsim_bp_Dict = dict([(key, value) for key, value in zip(semsim_bp['index'], semsim_bp['ss'])])

    # Process sequence similarity file
    seqsim['Protein_1'] = seqsim['Protein_1'].map(IDMappingDict)
    seqsim['Protein_2'] = seqsim['Protein_2'].map(IDMappingDict)
    seqsim['Protein_1'], seqsim['Protein_2'] = seqsim[['Protein_1', 'Protein_2']].min(axis=1), seqsim[['Protein_1', 'Protein_2']].max(axis=1)
    seqsim['index'] = seqsim[['Protein_1', 'Protein_2']].agg('_'.join, axis=1)
    identitiesDict = dict([(key, value) for key, value in zip(seqsim['index'], seqsim['Identities'])])
    positivesDict = dict([(key, value) for key, value in zip(seqsim['index'], seqsim['Positives'])])

    # # Find evalue cutoff
    # cutoff = distEvalue_paralogs['Evalue'].nlargest(2).iloc[1]
    # print(cutoff)

    cutoff = 1.0

    # Initialization of output
    result_all = []
    result_paralogs = []
    result_nonparalogs = []

    # Generate table
    for hub in hubProteins: 
        partnerList = get_partner_list(strInt, hub)
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        for i in range(len(partnerList)-1): 
            for j in range(i, len(partnerList)): 
                if partnerList[i] == partnerList[j]: 
                    continue
                index = '_'.join([min(partnerList[i], partnerList[j]), max(partnerList[i], partnerList[j])])
                degree_1, degree_2 = calculate_degree_of_interface_overlap(strInt_temp, hub, partnerList[i], partnerList[j])
                coexpressionScore = calculate_coexpression(coexpressionFile, partnerList[i], partnerList[j])
                if coexpressionScore == False: 
                    continue
                GIPS = calculate_GIPS(GI_profile, partnerList[i], partnerList[j])
                if GIPS == False: 
                    continue
                GI_share = calculate_GI_share(SGA_matrix, partnerList[i], partnerList[j])
                if GI_share == False: 
                    continue
                SS_mf = calculate_semantic_similarity(semsim_mf_Dict, partnerList[i], partnerList[j])
                if SS_mf == False: 
                    continue
                SS_bp = calculate_semantic_similarity(semsim_bp_Dict, partnerList[i], partnerList[j])
                if SS_bp == False: 
                    continue
                identities = calculate_sequence_similarity(identitiesDict, partnerList[i], partnerList[j])
                if identities == False: 
                    continue
                positives = calculate_sequence_similarity(positivesDict, partnerList[i], partnerList[j])
                if positives == False: 
                    continue
                evalue = calculate_evalue(distEvalue_all, partnerList[i], partnerList[j])
                if evalue == False: 
                    continue
                result_all.append((hub, partnerList[i], partnerList[j], degree_1, degree_2, coexpressionScore, GIPS, GI_share, SS_mf, SS_bp, identities, positives, evalue))
                if evalue <= cutoff: # Likely paralogous
                    result_paralogs.append((hub, partnerList[i], partnerList[j], degree_1, degree_2, coexpressionScore, GIPS, GI_share, SS_mf, SS_bp, identities, positives, evalue))
                else: 
                    result_nonparalogs.append((hub, partnerList[i], partnerList[j], degree_1, degree_2, coexpressionScore, GIPS, GI_share, SS_mf, SS_bp, identities, positives, evalue))
        # print("{} of {} hub proteins have been processed...".format(count, len(hubProteins)))
    print(result_all)
    print(result_paralogs)
    print(result_nonparalogs)
    result_all_df = pd.DataFrame(result_all, columns=['hub', 'partner_1', 'partner_2', 'degree_1', 'degree_2', 'coexpression', 'GIPS', 'SharedGIs', 'SemSim_mf', 'SemSim_bp', 'Identities', 'Positives', 'Evalue'])
    result_paralogs_df = pd.DataFrame(result_paralogs, columns=['hub', 'partner_1', 'partner_2', 'degree_1', 'degree_2', 'coexpression', 'GIPS', 'SharedGIs', 'SemSim_mf', 'SemSim_bp', 'Identities', 'Positives', 'Evalue'])
    result_nonparalogs_df = pd.DataFrame(result_nonparalogs, columns=['hub', 'partner_1', 'partner_2', 'degree_1', 'degree_2', 'coexpression', 'GIPS', 'SharedGIs', 'SemSim_mf', 'SemSim_bp', 'Identities', 'Positives', 'Evalue'])
    result_all_df.to_csv(OutPathAll, sep='\t', index=False)
    result_paralogs_df.to_csv(OutPathParalogs, sep='\t', index=False)
    result_nonparalogs_df.to_csv(OutPathNonparalogs, sep='\t', index=False)