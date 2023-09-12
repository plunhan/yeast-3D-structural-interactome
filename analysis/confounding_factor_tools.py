import pandas as pd 
import numpy as np
import pickle
import subprocess
import re
from scipy.stats import pearsonr, ranksums, gmean
from import_tools import (map_protein_id_to_locus_id, 
                          import_structural_interactome, 
                          import_structural_interactome_with_interacting_residues, 
                          import_hub_proteins, 
                          import_sequence_similarity, 
                          import_GI_profile, 
                          import_coexpression_data)
from interface_tools import (get_partner_list, 
                             calculate_degree_of_interface_overlap, 
                             calculate_GIPS)
from sequence_similarity_tools import (write_sequence_file, 
                                       retrieve_information, 
                                       generate_interactome_sequence_dict, 
                                       get_extended_interacting_residues_for_pair, 
                                       find_best_alignment)
from structural_similarity_tools import (calculate_structural_similarity, 
                                         find_best_TM_score)

# Essentiality 

def assign_essentiality(gene, essentialGeneList): 
    '''
    
    This function checks if a gene is essential or not. 

    Input: 

        gene (str): gene ordered locus name
        essentialGeneList (list): list of essential genes

    Return: 

        str: 'essential' for essential gene, 'non-essential' for non-essential gene
    '''
    if gene in essentialGeneList: 
        return 'essential'
    else: 
        return 'non-essential' 

def check_essentiality(InPathStrInt, 
                       InPathDegree, 
                       InPathIDMappingFile, 
                       InPathGIPS, 
                       InPathEssentialGeneList, 
                       OutPath): 
    '''
    
    This function checks if a gene is essential or not. 

    Input: 

        gene (str): gene ordered locus name
        essentialGeneList (list): list of essential genes

    Return: 

        str: 'essential' for essential gene, 'non-essential' for non-essential gene
    '''
    strInt = import_structural_interactome(InPathStrInt, InPathIDMappingFile)
    hubProteins, hubDict = import_hub_proteins(InPathDegree, InPathIDMappingFile)
    essentialGeneList = pickle.load(open(InPathEssentialGeneList, 'rb'))
    GI_profile = import_GI_profile(InPathGIPS)

    result_all = []
    for hub in hubProteins: 
        eHub = assign_essentiality(hub, essentialGeneList)
        partnerList = get_partner_list(strInt, hub)
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        for i in range(len(partnerList)-1): 
            for j in range(i, len(partnerList)): 
                if partnerList[i] == partnerList[j]: 
                    continue
                eP1 = assign_essentiality(partnerList[i], essentialGeneList)
                eP2 = assign_essentiality(partnerList[j], essentialGeneList)
                index = '_'.join([min(partnerList[i], partnerList[j]), max(partnerList[i], partnerList[j])])
                degree_1, degree_2 = calculate_degree_of_interface_overlap(strInt_temp, hub, partnerList[i], partnerList[j])
                GIPS = calculate_GIPS(GI_profile, partnerList[i], partnerList[j])
                if GIPS == False: 
                    continue
                result_all.append((hub, partnerList[i], partnerList[j], degree_1, degree_2, GIPS, eHub, eP1, eP2))
    result_all_df = pd.DataFrame(result_all, columns=['hub', 'partner_1', 'partner_2', 'degree_1', 'degree_2', 'GIPS', 'e-hub', 'e-partner_1', 'e-partner_2'])
    result_all_df.dropna(inplace=True)
    result_all_df.to_csv(OutPath, sep='\t', index=False, na_rep='NA')

###############################
###############################
###############################
##### Sequence similarity #####
###############################
###############################
###############################

## part I: whole gene sequence similarity

def generate_structural_interactome_protein_sequences(InPathProteome, InPathStrInt, OutPath): 
    '''
    
    This function returns a fasta file that contains sequences of proteins within the structural interactome. 

    Input: 

        InPathProteome (Path): path to fasta protein sequence file for all proteins in proteome
        InPathStrInt (Path): path to structural interactome file

    Output: 

        fasta: sequences of proteins within the structural interactome

    '''
    strInt = pd.read_table(InPathStrInt)
    strIntProteins = set(strInt[['Protein_1', 'Protein_2']].values.flatten().tolist())
    with open(OutPath, 'w') as o, open(InPathProteome, 'r') as i: 
        for line in i.readlines(): 
            if line[0] == '>': 
                protein = line.split('|')[1]
                if protein in strIntProteins: 
                    o.write(line)
            else: 
                if protein in strIntProteins: 
                    o.write(line)

def blast_interactome_protein_pairs(InPathIntSeq, 
                                    InPathStrInt, 
                                    InPathDegree, 
                                    InPathIDMappingFile, 
                                    InPathGIPS, 
                                    OutPathTempFile, 
                                    OutPathSeqQuery, 
                                    OutPathSeqSubject, 
                                    OutPath): 
    '''
    
    This function performs all-against-all blast for proteins within the structural interactome. 

    Input: 

        InPathIntSeq (Path): path to interactome sequences
        InPathStrInt (Path): path to structural interactome
        InPathDegree (Path): path to hub degree
        InPathIDMappingFile (Path): path to ID mapping file
        InPathGIPS (Path): path to genetic interaction matrix

    Output: 

        OutPathTempFile (Path): temporary file to save blast result
        OutPathSeqQuery (Path): temporary file to save query protein sequence
        OutPathSeqSubject (Path): temporary file to save subject protein sequence
        OutPath (Path): DataFrame - interactor pairs with hub protein, identities, positives, and e-value

    '''
    seqDict = generate_interactome_sequence_dict(InPathIntSeq)
    strInt = pd.read_table(InPathStrInt)
    degreeStrInt = pd.read_table(InPathDegree)
    hubProteins = list(set(degreeStrInt[degreeStrInt['Degree'] >= 2]['Protein']))
    IDMappingDict = map_protein_id_to_locus_id(InPathIDMappingFile)
    GI_profile = import_GI_profile(InPathGIPS)

    allPairList = []
    for hub in hubProteins: 
        partnerList = get_partner_list(strInt, hub)
        pairList = [(protein_1, protein_2) for idx, protein_1 in enumerate(partnerList) for protein_2 in partnerList[idx+1:]]
        allPairList += pairList
    allPartnerList = list(set(allPairList))

    seqSimDict = {}
    i = 1
    for pair in allPairList: 
        query, subject = pair[0], pair[1]
        write_sequence_file(query, subject, seqDict, OutPathSeqQuery, OutPathSeqSubject)
        cmd = ['blastp', '-query', str(OutPathSeqQuery), '-subject', str(OutPathSeqSubject), '-outfmt', '0', '-out', str(OutPathTempFile)]
        subprocess.run(cmd, capture_output = True)
        identities, positives, evalue = retrieve_information(OutPathTempFile)
        protein_1, protein_2 = IDMappingDict[query], IDMappingDict[subject]
        protein_1, protein_2 = min(protein_1, protein_2), max(protein_1, protein_2)
        index = '_'.join([protein_1, protein_2])
        seqSimDict[index] = (identities, positives, evalue)
        print('\r{} of {} ({:.2f}%%) protein pairs were processed...'.format(i, len(allPairList), i / len(allPairList) * 100), end = '')
        i += 1

    hubProteins = [IDMappingDict[protein] for protein in hubProteins]
    outputList = []
    strInt = import_structural_interactome(InPathStrInt, InPathIDMappingFile)
    count = 0
    for hub in hubProteins: 
        partnerList = get_partner_list(strInt, hub)
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        for i in range(len(partnerList)-1): 
            for j in range(i, len(partnerList)): 
                index = '_'.join([min(partnerList[i], partnerList[j]), max(partnerList[i], partnerList[j])])
                degree_1, degree_2 = calculate_degree_of_interface_overlap(strInt_temp, hub, partnerList[i], partnerList[j])
                GIPS = calculate_GIPS(GI_profile, partnerList[i], partnerList[j])
                if (GIPS == False) or (index not in seqSimDict.keys()): 
                    continue
                identities, positives, evalue = seqSimDict[index]
                outputList.append((hub, partnerList[i], partnerList[j], degree_1, degree_2, GIPS, identities, positives, evalue))
        print('\r{} of {} ({:.2f}%%) hubs were processed...'.format(count, len(hubProteins), count / len(hubProteins) * 100), end = '')

    outputDataFrame = pd.DataFrame(outputList, columns=['hub', 'Protein_1', 'Protein_2', 'ORC', 'ORR', 'GIPS', 'Identities', 'Positives', 'Evalue'])
    outputDataFrame.to_csv(OutPath, sep='\t', na_rep = 'NA')

## part II: sequence similarity of interfacial residues

def produce_sequence_similarity_covered_interfaces(InPathStrInt, 
                                                   InPathDegree, 
                                                   InPathIDMappingFile, 
                                                   InPathGIPS, 
                                                   InPathProteome, 
                                                   hubCutOff, 
                                                   TempFile_1, 
                                                   TempFile_2, 
                                                   blastTempFile, 
                                                   alignment_cutoff,  
                                                   OutPath): 
    '''
    
    This function maps sequence similarity (BLAST E-value) to each pair of interactors on the same target protein 
    within the structural interactome. 

    Input: 

        InPathStrInt (Path): path to structural interactome
        InPathDegree (Path): path to hub degree
        InPathIDMappingFile (Path): path to ID mapping file
        InPathGIPS (Path): path to genetic interaction matrix
        InPathProteome (Path): path to protein sequence file
        hubCutOff (int): degree of PPI which defines whether a protein is a target protein
        TempFile_1 (Path): path to temporarily save the sequence of the first protein during the BLAST pair-wise
                           sequence alignment
        TempFile_2 (Path): path to the file that temporarily saves the sequence of the second protein during the BLAST pair-wise
                           sequence alignment
        blastTempFile (Path): path to the file that temporarily saves the output of the BLAST pair-wise sequence alignment
        alignmnet_cutoff (int): the fraction cutoff of interfacial residues on each protein that are aligned by BLAST

    Output: 

        OutPath (Path): DataFrame - interactor pairs with target protein, interfacial overlap, genetic interaction 
                        profile similarity, and E-value

    '''
    strInt = pickle.load(open(InPathStrInt, 'rb'))
    degreeStrInt = pd.read_table(InPathDegree)
    hubProteins = list(set(degreeStrInt[degreeStrInt['Degree'] >= hubCutOff]['Protein']))
    IDMappingDict = map_protein_id_to_locus_id(InPathIDMappingFile)
    GI_profile = import_GI_profile(InPathGIPS)
    allProtSeq = pickle.load(open(InPathProteome, 'rb'))

    resultList = []
    count_hub = 0
    for hub in hubProteins: 
        count_hub += 1
        partnerList = get_partner_list(strInt, hub)
        partnerPairList = [(partner_1, partner_2) for idx, partner_1 in enumerate(partnerList) for partner_2 in partnerList[idx+1:]]
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        count_pairs = 0
        for partnerPair in partnerPairList: 
            count_pairs += 1
            protein_1, protein_2 = partnerPair
            if protein_1 == protein_2: 
                continue
            interactingRes1, interactingRes2, ORC, ORCR = get_extended_interacting_residues_for_pair(hub, 
                                                                                                     protein_1, 
                                                                                                     protein_2, 
                                                                                                     strInt_temp)
            evalue = find_best_alignment(hub, 
                                         protein_1, 
                                         protein_2, 
                                         interactingRes1, 
                                         interactingRes2, 
                                         allProtSeq, 
                                         TempFile_1, 
                                         TempFile_2, 
                                         blastTempFile, 
                                         alignment_cutoff, 
                                         count_hub, 
                                         len(hubProteins), 
                                         count_pairs, 
                                         len(partnerPairList))
            protein_1, protein_2 = IDMappingDict[protein_1], IDMappingDict[protein_2]
            GIPS = calculate_GIPS(GI_profile, protein_1, protein_2)
            if not GIPS: 
                continue
            resultList.append( (IDMappingDict[hub], protein_1, protein_2, ORC, ORCR, GIPS, evalue) )
    resultDataFrame = pd.DataFrame(resultList, columns = ['hub', 'Protein_1', 'Protein_2', 'ORC', 'ORCR', 'GIPS', 'E-value'])
    resultDataFrame.to_csv(OutPath, sep = '\t', na_rep = 'NA')

#################################
#################################
#################################
##### Structural similarity #####
#################################
#################################
#################################

# structural similarity of whole proteins

def produce_structural_similarity_whole_protein(InPathStrInt, 
                                                InPathDegree, 
                                                InPathIDMappingFile, 
                                                InPathGIPS, 
                                                ppiModelDir, 
                                                hubCutOff, 
                                                TempFile_1, 
                                                TempFile_2, 
                                                TMalign_file, 
                                                TMalign_path, 
                                                OutPath): 
    '''
    
    This function maps structural similarity (TM-score generated by TM-align) to each pair of interactors on 
    the same target protein within the structural interactome. 

    Input: 

        InPathStrInt (Path): path to structural interactome
        InPathDegree (Path): path to hub degree
        InPathIDMappingFile (Path): path to ID mapping file
        InPathGIPS (Path): path to genetic interaction matrix
        ppiModelDir (Path): path to directory of PPI models
        hubCutOff (int): degree of PPI which defines whether a protein is a target protein
        TempFile_1 (Path): path to the file that temporarily saves the sequence of the first protein during the TM-align
                           structural alignment
        TempFile_2 (Path): path to the file that temporarily saves the sequence of the second protein during the TM-align
                           structural alignment
        TMalign_file (Path): path to the file that temporarily saves the output of the TM-align structural alignment
        TMalign_path (Path): path to TM-align software

    Output: 

        OutPath (Path): DataFrame - interactor pairs with target protein, interfacial overlap, genetic interaction 
                        profile similarity, and TM-score

    '''
    strInt = pd.read_table(InPathStrInt)
    degreeStrInt = pd.read_table(InPathDegree)
    hubProteins = list(set(degreeStrInt[degreeStrInt['Degree'] >= hubCutOff]['Protein']))
    IDMappingDict = map_protein_id_to_locus_id(InPathIDMappingFile)
    GI_profile = import_GI_profile(InPathGIPS)

    resultList = []
    count_hub = 0
    for hub in hubProteins: 
        count_hub += 1
        partnerList = get_partner_list(strInt, hub)
        partnerPairList = [(partner_1, partner_2) for idx, partner_1 in enumerate(partnerList) for partner_2 in partnerList[idx+1:]]
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        count_pairs = 0
        for partnerPair in partnerPairList: 
            count_pairs += 1
            protein_1, protein_2 = partnerPair
            if protein_1 == protein_2: 
                continue
            strSim = calculate_structural_similarity(hub, 
                                                     protein_1, 
                                                     protein_2, 
                                                     TempFile_1, 
                                                     TempFile_2, 
                                                     TMalign_file, 
                                                     TMalign_path, 
                                                     ppiModelDir, 
                                                     count_hub, 
                                                     len(hubProteins),
                                                     count_pairs, 
                                                     len(partnerPairList))
            ORC = calculate_degree_of_interface_overlap(strInt_temp, hub, protein_1, protein_2)
            protein_1, protein_2 = IDMappingDict[protein_1], IDMappingDict[protein_2]
            GIPS = calculate_GIPS(GI_profile, protein_1, protein_2)
            if not GIPS: 
                continue
            resultList.append( (IDMappingDict[hub], protein_1, protein_2, ORC, GIPS, strSim) )
    resultDataFrame = pd.DataFrame(resultList, columns = ['hub', 'Protein_1', 'Protein_2', 'ORC', 'GIPS', 'TM-score'])
    resultDataFrame.to_csv(OutPath, sep = '\t', na_rep = 'NA')

# structual similarity of interface

def produce_structural_similarity_extended_interfaces(InPathStrInt, 
                                                      InPathDegree, 
                                                      InPathIDMappingFile, 
                                                      InPathGIPS, 
                                                      InPathProteome, 
                                                      ppiModelDir, 
                                                      hubCutOff, 
                                                      TempFile_1, 
                                                      TempFile_2, 
                                                      TMalign_file, 
                                                      TMalign_path, 
                                                      strSim_cutoff, 
                                                      OutPath): 
    '''
    This function maps structural similarity (TM-score generated by TM-align) to each pair of interactors on the same target protein 
    within the structural interactome. 

    Input: 

        InPathStrInt (Path): path to structural interactome
        InPathDegree (Path): path to hub degree
        InPathIDMappingFile (Path): path to ID mapping file
        InPathGIPS (Path): path to genetic interaction matrix
        InPathProteome (Path): path to protein sequence file
        ppiModelDir (Path): path to directory of PPI models
        hubCutOff (int): degree of PPI which defines whether a protein is a target protein
        TempFile_1 (Path): path to the file that temporarily saves the sequence of the first protein during the TM-align
                           structural alignment
        TempFile_2 (Path): path to the file that temporarily saves the sequence of the second protein during the TM-align
                           structural alignment
        TMalign_file (Path): path to the file that temporarily saves the output of the TM-align structural alignment
        TMalign_path (Path): path to TM-align software
        strSim_cutoff (int): the fraction cutoff of interfacial residues on each protein that are aligned by TM-align

    Output: 

        OutPath (Path): DataFrame - interactor pairs with target protein, interfacial overlap, genetic interaction 
                        profile similarity, and TM-scores with 10%, 30%, and 50% of interfacial sequence segment aligned

    '''
    strInt = pickle.load(open(InPathStrInt, 'rb'))
    degreeStrInt = pd.read_table(InPathDegree)
    hubProteins = list(set(degreeStrInt[degreeStrInt['Degree'] >= hubCutOff]['Protein']))
    IDMappingDict = map_protein_id_to_locus_id(InPathIDMappingFile)
    GI_profile = import_GI_profile(InPathGIPS)
    allProtSeq = pickle.load(open(InPathProteome, 'rb'))

    resultList = []
    count_hub = 0
    for hub in hubProteins: 
        count_hub += 1
        partnerList = get_partner_list(strInt, hub)
        partnerPairList = [(partner_1, partner_2) for idx, partner_1 in enumerate(partnerList) for partner_2 in partnerList[idx+1:]]
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        count_pairs = 0
        for partnerPair in partnerPairList: 
            count_pairs += 1
            protein_1, protein_2 = partnerPair
            if protein_1 == protein_2: 
                continue
            interactingRes1, interactingRes2, ORC, ORCR = get_extended_interacting_residues_for_pair(hub, 
                                                                                                     protein_1, 
                                                                                                     protein_2, 
                                                                                                     strInt_temp)
            TMscore_50, TMscore_30, TMscore_10 = find_best_TM_score(hub, 
                                                                    protein_1, 
                                                                    protein_2, 
                                                                    interactingRes1, 
                                                                    interactingRes2, 
                                                                    allProtSeq, 
                                                                    ppiModelDir, 
                                                                    TempFile_1, 
                                                                    TempFile_2, 
                                                                    TMalign_file, 
                                                                    TMalign_path, 
                                                                    count_hub, 
                                                                    len(hubProteins), 
                                                                    count_pairs, 
                                                                    len(partnerPairList), 
                                                                    strSim_cutoff)
            protein_1, protein_2 = IDMappingDict[protein_1], IDMappingDict[protein_2]
            GIPS = calculate_GIPS(GI_profile, protein_1, protein_2)
            if not GIPS: 
                continue
            resultList.append( (IDMappingDict[hub], protein_1, protein_2, ORC, ORCR, GIPS, TMscore_50, TMscore_30, TMscore_10) )
    resultDataFrame = pd.DataFrame(resultList, columns = ['hub', 'Protein_1', 'Protein_2', 'ORC', 'ORCR', 'GIPS', 'TM-score_50', 'TM-score_30', 'TM-score_10'])
    resultDataFrame.to_csv(OutPath, sep = '\t', na_rep = 'NA')