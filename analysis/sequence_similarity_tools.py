#----------------------------------------------------------------------------------------
# Modules for computations on sequence similarity of two proteins
#----------------------------------------------------------------------------------------

import pickle
import itertools
import pandas as pd
import subprocess
from import_tools import map_protein_id_to_locus_id


def delete_file(file): 
    if file.is_file(): 
        subprocess.call(['rm', str(file)])

# used for whole gene sequence similarity
def generate_interactome_sequence_dict(InPath):
    seqDict = {}
    curProt = ''
    curSeq = ''
    with open(InPath, 'r') as f: 
        for line in f: 
            if line.startswith('>') and curProt == '': 
                curProt = line.rstrip().split('>')[1]
            elif line.startswith('>') and curProt != '': 
                seqDict[curProt] = curSeq
                curProt = line.rstrip().split('>')[1]
                curSeq = ''
            else: 
                curSeq = line.rstrip()
        seqDict[curProt] = curSeq
    return seqDict

def retrieve_information(Path): # identities, positives, and e-value
    with open(Path, 'r') as f: 
        string = f.read()
        if 'No hits found' in string: 
            return 0, 0, 1000
        position_identities = string.find('Identities') + len('Identities') + 3
        identities = float(string[position_identities:].split('/')[0]) / float(string[position_identities:].split('/')[1].split(' ')[0])
        position_positives = string.find('Positives') + len('Positives') + 3
        positives = float(string[position_positives:].split('/')[0]) / float(string[position_positives:].split('/')[1].split(' ')[0])
        evalue = float(string[string.find('Expect'):].split(',')[0].split(' ')[-1])
    return identities, positives, evalue

def write_sequence_file(query, 
                        subject, 
                        seqDict, 
                        OutPathSeqQuery, 
                        OutPathSeqSubject): 
    seqQuery = seqDict[query]
    with open(OutPathSeqQuery, 'w') as f: 
        f.write('>' + query + '\n' + seqQuery)
    seqSubject = seqDict[subject]
    with open(OutPathSeqSubject, 'w') as f: 
        f.write('>' + subject + '\n' + seqSubject)

# used for mapping Pfam domains
def get_pfam_list(protein, Pfam): 
    '''
    
    This function returns a list of Pfam names. 

    Input: 

        protein (str): gene ordered locus name
        Pfam (DataFrame): dataframe of Pfam mapping result

    Return: 

        str: 'essential' for essential gene, 'non-essential' for non-essential gene
    '''
    return list(Pfam[Pfam['seqid'] == protein]['hmm.acc'])

def get_interacting_residues_for_pair(hub, 
                                      interactor1, 
                                      interactor2, 
                                      strInt): 
    '''

    Input: 

        hub (str): A specific hub protein. 
        interactor1 (str): One interactor of the hub protein. 
        interactor2 (str): Another interactor of the hub protein. 
        strInt (DataFrame): Structurally-resolved PPIs for the hub protein. 

    Return: 

        List, list - Residues interacting with "overlapping residues on hub protein". 

    '''
    ORC, ORCR = get_overlapping_residues(hub, 
                                         interactor1, 
                                         interactor2, 
                                         strInt)
    interactingRes1 = get_interacting_residues_for_interactor(hub, 
                                                              interactor1, 
                                                              strInt)
    interactingRes2 = get_interacting_residues_for_interactor(hub, 
                                                              interactor2, 
                                                              strInt)
    return interactingRes1, interactingRes2, ORC, ORCR

def get_extended_interacting_residues_for_pair(hub, 
                                               interactor1, 
                                               interactor2, 
                                               strInt): 
    '''

    Input: 

        hub (str): A specific hub protein. 
        interactor1 (str): One interactor of the hub protein. 
        interactor2 (str): Another interactor of the hub protein. 
        strInt (DataFrame): Structurally-resolved PPIs for the hub protein. 

    Return: 

        List, list - Residues interacting with "overlapping residues on hub protein". 

    '''
    ORC, ORCR = get_overlapping_residues(hub, 
                                         interactor1, 
                                         interactor2, 
                                         strInt)
    interactingRes1 = get_extended_interacting_residues_for_interactor(strInt, 
                                                                       interactor1, 
                                                                       strInt)
    interactingRes2 = get_extended_interacting_residues_for_interactor(strInt, 
                                                                       interactor2, 
                                                                       strInt)
    return interactingRes1, interactingRes2, ORC, ORCR

def get_overlapping_residues(hub, 
                             interactor1, 
                             interactor2, 
                             strInt): 
    '''

    Input: 

        hub (str): A specific hub protein. 
        interactor1 (str): One interactor of the hub protein. 
        interactor2 (str): Another interactor of the hub protein. 
        strInt (DataFrame): Structurally-resolved PPIs for the hub protein. 

    Return: 

        List - Overlapping residues on hub protein. 

    '''

    strInt1 = strInt[(strInt.Protein_1 == interactor1) 
                   | (strInt.Protein_2 == interactor1)].iloc[0]

    if strInt1.Protein_1 == interactor1: 
        residues1 = strInt1['Interfaces'].split('+')[1].split(',')
    elif strInt1.Protein_2 == interactor1: 
        residues1 = strInt1['Interfaces'].split('+')[0].split(',')
    residues1 = [int(residue) for residue in residues1]

    strInt2 = strInt[(strInt.Protein_1 == interactor2) 
                   | (strInt.Protein_2 == interactor2)].iloc[0]

    if strInt2.Protein_1 == interactor2: 
        residues2 = strInt2['Interfaces'].split('+')[1].split(',')
    elif strInt2.Protein_2 == interactor2: 
        residues2 = strInt2['Interfaces'].split('+')[0].split(',')
    residues2 = [int(residue) for residue in residues2]

    OR = list(set(residues1).intersection(set(residues2)))
    allOR = list(set(residues1).union(set(residues2)))

    return len(OR), len(OR)/len(allOR)

def get_interacting_residues_for_interactor(strInt_temp, 
                                            interactor, 
                                            strInt): 
    '''

    Input: 

        hub (str): A specific hub protein. 
        interactor1 (str): One interactor of the hub protein. 
        interactor2 (str): Another interactor of the hub protein. 
        strInt (DataFrame): Structurally-resolved PPIs for the hub protein. 

    Return: 

        List - Residues interacting with "overlapping residues on hub protein". 

    '''
    strInt_temp = strInt[(strInt.Protein_1 == interactor)
                       | (strInt.Protein_2 == interactor)].iloc[0]
    overlappingResIndex = 1 if strInt_temp.Protein_1 == interactor else 0
    if overlappingResIndex: 
        return list(set([intRes for intRes, overRes in strInt_temp['Interacting residue pairs']]))
    else: 
        return list(set([intRes for overRes, intRes in strInt_temp['Interacting residue pairs']]))

def get_extended_interacting_residues_for_interactor(strInt_temp, 
                                                     interactor, 
                                                     strInt): 
    '''

    Input: 

        hub (str): A specific hub protein. 
        interactor1 (str): One interactor of the hub protein. 
        interactor2 (str): Another interactor of the hub protein. 
        strInt (DataFrame): Structurally-resolved PPIs for the hub protein. 

    Return: 

        List - Sequence of the smallest protein fragments containing all interfacial residues. 

    '''
    strInt_temp = strInt[(strInt.Protein_1 == interactor)
                       | (strInt.Protein_2 == interactor)].iloc[0]
    overlappingResIndex = 1 if strInt_temp.Protein_1 == interactor else 0
    if overlappingResIndex: 
        res = list(set([intRes for intRes, overRes in strInt_temp['Interacting residue pairs']]))
    else: 
        res = list(set([intRes for overRes, intRes in strInt_temp['Interacting residue pairs']]))
    return list(range(min(res), max(res)+1))

def get_competing_residues_on_Pfam(interactor1, 
                                   interactor2, 
                                   interactingRes1, 
                                   interactingRes2, 
                                   sharedPfams, 
                                   Pfam): 
    '''

    Input: 
        
        interactor1 (str): Name of interactor 1. 
        interactor2 (str): Name of interactor 2. 
        interactingRes1 (list): A list of interactor 1 residues interacting with overlapping residues. 
        interactingRes2 (list): A list of interactor 2 residues interacting with overlapping residues. 
        sharedPfams (list): A list of Pfam families shared by interactors 1 and 2. 
        Pfam (DataFrame): Pfam scan results. 

    Return: 

        float, float: The percentage of interacting residues belonging to shared Pfam families. 

    '''

    Pfam1 = Pfam[Pfam.seqid == interactor1]
    Pfam2 = Pfam[Pfam.seqid == interactor2]
    PfamDict1 = generate_Pfam_dict_for_interactor(Pfam1, interactor1)
    PfamDict2 = generate_Pfam_dict_for_interactor(Pfam2, interactor2)
    PfamDict1IntRes = dict([(key, value) for key, value in PfamDict1.items() if key in interactingRes1])
    PfamDict2IntRes = dict([(key, value) for key, value in PfamDict2.items() if key in interactingRes2])
    Pfam1IntRes = set().union(*PfamDict1IntRes.values())
    Pfam2IntRes = set().union(*PfamDict2IntRes.values())
    realSharedPfams = Pfam1IntRes.intersection(Pfam2IntRes)
    result1 = [res for res, Pfam in PfamDict1.items() if (res in interactingRes1) and len(Pfam.intersection(realSharedPfams))]
    result2 = [res for res, Pfam in PfamDict2.items() if (res in interactingRes2) and len(Pfam.intersection(realSharedPfams))]
    return len(result1) / len(interactingRes1), len(result2) / len(interactingRes2)

def generate_Pfam_dict_for_interactor(Pfam, interactor): 
    '''

    Input: 

        Pfam (DataFrame): Pfam family mapping file for yeast. 
        interactor (str): Target interactor. 

    Return: 

        dict: (key, value) - (residue, a set of Pfam families map to the residue)
    '''
    PfamDict = {}
    for i, row in Pfam.iterrows():
        if row['seqid'] != interactor: 
            continue
        PfamID = row['hmm.acc']
        residues = list(range(row['alignment.start'], row['alignment.end']+1))
        for residue in residues: 
            if residue in PfamDict.keys():
                PfamDict[residue].add(PfamID)
            else: 
                PfamDict[residue] = set([PfamID])
    return PfamDict

def calculate_local_sequence_similarity(hub, 
                                        protein_1, 
                                        protein_2, 
                                        PfamP1, 
                                        PfamP2, 
                                        PfamToResDict, 
                                        allProtSeq, 
                                        TempFile_1, 
                                        TempFile_2, 
                                        OutPathTempFile, 
                                        count_hub, 
                                        total_hub, 
                                        count_pairs, 
                                        total_pairs): 
    allSeq_1 = allProtSeq[protein_1]
    allSeq_2 = allProtSeq[protein_2]
    res_1 = []
    res_2 = []
    for Pfam in PfamP1: 
        res_1 += PfamToResDict[protein_1][Pfam]
    for Pfam in PfamP2: 
        res_2 += PfamToResDict[protein_2][Pfam]
    res_1 = set(list(range(min(res_1), max(res_1)+1)))
    res_2 = set(list(range(min(res_2), max(res_2)+1)))
    seq_1 = ''.join([res for idx, res in enumerate(allSeq_1) if idx+1 in res_1])
    print(seq_1)
    seq_2 = ''.join([res for idx, res in enumerate(allSeq_2) if idx+1 in res_2])
    print(seq_2)
    with open(TempFile_1, 'w') as f: 
        f.write('>' + protein_1 + '\n' + seq_1)
    with open(TempFile_2, 'w') as f: 
        f.write('>' + protein_2 + '\n' + seq_2)
    cmd = ['blastp', '-query', str(TempFile_1), '-subject', str(TempFile_2), '-outfmt', '0', '-out', str(OutPathTempFile)]
    subprocess.run(cmd, capture_output = True)
    identities, positives, evalue = retrieve_information(OutPathTempFile)
    delete_file(TempFile_1)
    delete_file(TempFile_2)
    delete_file(OutPathTempFile)
    # print('Hub protein: {}'.format(hub))
    # print('    {} of {} hub proteins have been processed. '.format(count_hub, total_hub))
    # print('    Calculating sequence similarity for protein pairs: ')
    # print('         {} + {}'.format(protein_1, protein_2))
    # print('         {} of {} interactor pairs have been processed. '.format(count_pairs, total_pairs))
    # print('    E-value for {} and {} is: '.format(protein_1, protein_2))
    # print('         {}'.format(evalue))
    return identities, positives, evalue

def identical_characters(string_1, string_2): 
    return sum(l1 == l2 for l1, l2 in zip(string_1, string_2))

def fraction_interface(blastSeq, intSeq): 
    return len([res for res in intSeq if res in blastSeq])

def find_best_alignment(hub, 
                        protein_1, 
                        protein_2, 
                        res_1, 
                        res_2, 
                        allProtSeq, 
                        TempFile_1, 
                        TempFile_2, 
                        OutPathTempFile, 
                        alignment_cutoff, 
                        count_hub, 
                        total_hub, 
                        count_pairs, 
                        total_pairs): 
    evalue_output = -1
    allSeq_1 = allProtSeq[protein_1]
    allSeq_2 = allProtSeq[protein_2]
    with open(TempFile_1, 'w') as f: 
        f.write('>' + protein_1 + '\n' + allSeq_1)
    with open(TempFile_2, 'w') as f: 
        f.write('>' + protein_2 + '\n' + allSeq_2)
    cmd = ['blastp', '-query', str(TempFile_1), '-subject', str(TempFile_2), '-outfmt', '6 qseqid sseqid qstart qend sstart send qseq sseq evalue', '-out', str(OutPathTempFile)]
    subprocess.run(cmd, capture_output = True)
    blast_result = pd.read_table(OutPathTempFile, header=None, names=['qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq', 'evalue'])
    if len(blast_result) == 0: 
        return 100
    blast_result = blast_result.sort_values(by=['evalue'])
    flag = 0
    for index, row in blast_result.iterrows(): 
        if flag: 
            break
        query_res = list(range(row['qstart'], row['qend']+1))
        subject_res = list(range(row['sstart'], row['send']+1))
        if fraction_interface(query_res, res_1) >= alignment_cutoff and fraction_interface(subject_res, res_2) >= alignment_cutoff: 
            flag = 1
            evalue_output = row['evalue']
    print('Hub protein: {}'.format(hub))
    print('    {} of {} hub proteins have been processed. '.format(count_hub, total_hub))
    print('    Calculating sequence similarity for protein pairs: ')
    print('         {} + {}'.format(protein_1, protein_2))
    print('         {} of {} interactor pairs have been processed. '.format(count_pairs, total_pairs))
    print('    E-value for {} and {} is: '.format(protein_1, protein_2))
    print('         {}'.format(evalue_output))
    return evalue_output