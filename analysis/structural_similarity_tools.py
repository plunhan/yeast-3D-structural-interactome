import subprocess
import re
from pathlib import Path
from sequence_similarity_tools import (get_interacting_residues_for_interactor, delete_file)
import math
import colorama
from TMalign_tools import check_mapping_status_for_interface

def calculate_structural_similarity(hub, 
                                    protein_1, 
                                    protein_2, 
                                    TempFile_1, 
                                    TempFile_2, 
                                    TMalign_file, 
                                    TMalign_path, 
                                    ppiModelDir, 
                                    count_hub, 
                                    total_hub, 
                                    count_pairs, 
                                    total_pairs): 
    file_1_1 = 'pdb' + hub + '=' + protein_1 + '.ent'
    path_1_1 = ppiModelDir / file_1_1
    file_1_2 = 'pdb' + protein_1 + '=' + hub + '.ent'
    path_1_2 = ppiModelDir / file_1_2
    file_2_1 = 'pdb' + hub + '=' + protein_2 + '.ent'
    path_2_1 = ppiModelDir / file_2_1
    file_2_2 = 'pdb' + protein_2 + '=' + hub + '.ent'
    path_2_2 = ppiModelDir / file_2_2
    if path_1_1.is_file(): 
        path_1 = path_1_1
        pos_1 = 'B'
    elif path_1_2.is_file(): 
        path_1 = path_1_2
        pos_1 = 'A'
    else: 
        return -1
    if path_2_1.is_file(): 
        path_2 = path_2_1
        pos_2 = 'B'
    elif path_2_2.is_file(): 
        path_2 = path_2_2
        pos_2 = 'A'
    else: 
        return -1
    TMscore = calculate_TM_score_whole_protein(path_1, 
                                               pos_1, 
                                               path_2, 
                                               pos_2, 
                                               TempFile_1, 
                                               TempFile_2, 
                                               TMalign_file, 
                                               TMalign_path, 
                                               ppiModelDir)
    print('Hub protein: {}'.format(hub))
    print('    {} of {} hub proteins have been processed. '.format(count_hub, total_hub))
    print('    Calculating structural similarity for protein pairs: ')
    print('         {} + {}'.format(protein_1, protein_2))
    print('         {} of {} interactor pairs have been processed. '.format(count_pairs, total_pairs))
    print('    TM-score for {} and {} is: '.format(protein_1, protein_2))
    print('         {}'.format(TMscore))
    return TMscore

def calculate_TM_score_whole_protein(path_1, 
                                     pos_1, 
                                     path_2, 
                                     pos_2, 
                                     TempFile_1, 
                                     TempFile_2, 
                                     TMalign_file, 
                                     TMalign_path, 
                                     ppiModelDir): 
    with open(path_1, 'r') as fin, open(TempFile_1, 'w') as fout: 
        for line in fin.readlines(): 
            if line.startswith('ATOM'): 
                if line[21] == pos_1: 
                    fout.write(line)
            elif line.startswith('TER'): 
                if line[21] == pos_1: 
                    fout.write(line)
            elif line.startswith('END'):
                fout.write(line)
    with open(path_2, 'r') as fin, open(TempFile_2, 'w') as fout: 
        for line in fin.readlines(): 
            if line.startswith('ATOM'): 
                if line[21] == pos_2: 
                    fout.write(line)
            elif line.startswith('TER'): 
                if line[21] == pos_2: 
                    fout.write(line)
            elif line.startswith('END'):
                fout.write(line)
    cmd = [str(TMalign_path), str(TempFile_1), str(TempFile_2), '-a', 'T']
    summary = subprocess.run(cmd, capture_output = True)
    TMscore = retrieve_TM_score(summary.stdout)
    delete_file(TMalign_file)
    delete_file(TempFile_1)
    delete_file(TempFile_2)
    return TMscore

def retrieve_TM_score(line): 
    try: 
        TMscore = float(re.findall(r"\b\d+\.\d+\s\(", str(line))[2].split(' ')[0])
        return TMscore
    except: 
        return -1

def produce_interacting_residues_Pfam_annotation(Pfam): 
    PfamDict = {}
    for index, row in Pfam.iterrows():
        if row['seqid'] not in PfamDict: 
            PfamDict[row['seqid']] = []
        PfamDict[row['seqid']] += [(res, row['hmm.acc']) for res in range(row['alignment.start'], row['alignment.end'] + 1)]
    for key, value in PfamDict.items(): 
        PfamDict[key] = dict(value)
    return PfamDict

def produce_Pfam_pos_interacting_residues(Pfam): 
    PfamDict = {}
    for index, row in Pfam.iterrows(): 
        if row['seqid'] not in PfamDict: 
            PfamDict[row['seqid']] = {}
        PfamDict[row['seqid']][row['hmm.acc']] = list(range(row['alignment.start'], row['alignment.end']+1))
    return PfamDict

def get_interfacial_pfam_list(resToPfamDict, 
                              protein, 
                              interactingRes): 
    if protein not in resToPfamDict.keys():
        return []
    PfamSet = set()
    for res in interactingRes: 
        if res in resToPfamDict[protein].keys(): 
            PfamSet.add(resToPfamDict[protein][res])
    return list(PfamSet)

def calculate_local_structural_similarity(hub, 
                                          protein_1, 
                                          protein_2, 
                                          Pfam_1, 
                                          Pfam_2, 
                                          PfamToResDict, 
                                          TempFile_1, 
                                          TempFile_2, 
                                          TMalign_file, 
                                          TMalign_path, 
                                          ppiModelDir, 
                                          count_hub, 
                                          total_hub, 
                                          count_pairs, 
                                          total_pairs): 
    file_1_1 = 'pdb' + hub + '=' + protein_1 + '.ent'
    path_1_1 = ppiModelDir / file_1_1
    file_1_2 = 'pdb' + protein_1 + '=' + hub + '.ent'
    path_1_2 = ppiModelDir / file_1_2
    file_2_1 = 'pdb' + hub + '=' + protein_2 + '.ent'
    path_2_1 = ppiModelDir / file_2_1
    file_2_2 = 'pdb' + protein_2 + '=' + hub + '.ent'
    path_2_2 = ppiModelDir / file_2_2
    if path_1_1.is_file(): 
        path_1 = path_1_1
        pos_1 = 'B'
    elif path_1_2.is_file(): 
        path_1 = path_1_2
        pos_1 = 'A'
    else: 
        return -1
    if path_2_1.is_file(): 
        path_2 = path_2_1
        pos_2 = 'B'
    elif path_2_2.is_file(): 
        path_2 = path_2_2
        pos_2 = 'A'
    else: 
        return -1
    TMscore = calculate_TM_score_local(path_1, 
                                       pos_1, 
                                       path_2, 
                                       pos_2, 
                                       protein_1, 
                                       protein_2, 
                                       Pfam_1, 
                                       Pfam_2, 
                                       PfamToResDict, 
                                       TempFile_1, 
                                       TempFile_2, 
                                       TMalign_file, 
                                       TMalign_path, 
                                       ppiModelDir)
    print('Hub protein: {}'.format(hub))
    print('    {} of {} hub proteins have been processed. '.format(count_hub, total_hub))
    print('    Calculating structural similarity for protein pairs: ')
    print('         {} + {}'.format(protein_1, protein_2))
    print('         {} of {} interactor pairs have been processed. '.format(count_pairs, total_pairs))
    print('    TM-score for {} and {} is: '.format(protein_1, protein_2))
    print('         {}'.format(TMscore))
    return TMscore

def calculate_TM_score_local(path_1, 
                             pos_1, 
                             path_2, 
                             pos_2, 
                             protein_1, 
                             protein_2, 
                             Pfam_1, 
                             Pfam_2, 
                             PfamToResDict, 
                             TempFile_1, 
                             TempFile_2, 
                             TMalign_file, 
                             TMalign_path, 
                             ppiModelDir): 
    res_1 = PfamToResDict[protein_1][Pfam_1]
    res_2 = PfamToResDict[protein_2][Pfam_2]
    with open(path_1, 'r') as fin, open(TempFile_1, 'w') as fout: 
        for line in fin.readlines(): 
            if line.startswith('ATOM') and line[21] == pos_1 and int(line[22:26]) in res_1: 
                fout.write(line)
    with open(path_2, 'r') as fin, open(TempFile_2, 'w') as fout: 
        for line in fin.readlines(): 
            if line.startswith('ATOM') and line[21] == pos_1 and int(line[22:26]) in res_2: 
                fout.write(line)
    cmd = [str(TMalign_path), str(TempFile_1), str(TempFile_2), '-a', 'T']
    summary = subprocess.run(cmd, capture_output = True)
    TMscore = retrieve_TM_score(summary.stdout)
    delete_file(TMalign_file)
    delete_file(TempFile_1)
    delete_file(TempFile_2)
    return TMscore

def find_best_TM_score(hub, 
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
                       total_hub, 
                       count_pairs, 
                       total_pairs, 
                       strSim_cutoff): 
    file_1_1 = 'pdb' + hub + '=' + protein_1 + '.ent'
    path_1_1 = ppiModelDir / file_1_1
    file_1_2 = 'pdb' + protein_1 + '=' + hub + '.ent'
    path_1_2 = ppiModelDir / file_1_2
    file_2_1 = 'pdb' + hub + '=' + protein_2 + '.ent'
    path_2_1 = ppiModelDir / file_2_1
    file_2_2 = 'pdb' + protein_2 + '=' + hub + '.ent'
    path_2_2 = ppiModelDir / file_2_2
    if path_1_1.is_file(): 
        path_1 = path_1_1
        pos_1 = 'B'
    elif path_1_2.is_file(): 
        path_1 = path_1_2
        pos_1 = 'A'
    else: 
        return -1
    if path_2_1.is_file(): 
        path_2 = path_2_1
        pos_2 = 'B'
    elif path_2_2.is_file(): 
        path_2 = path_2_2
        pos_2 = 'A'
    else: 
        return -1
    TMscore_50, TMscore_30, TMscore_10 = calculate_TM_score_local_extended(path_1, 
                                                                           pos_1, 
                                                                           path_2, 
                                                                           pos_2, 
                                                                           protein_1, 
                                                                           protein_2, 
                                                                           allProtSeq, 
                                                                           interactingRes1, 
                                                                           interactingRes2, 
                                                                           TempFile_1, 
                                                                           TempFile_2, 
                                                                           TMalign_file, 
                                                                           TMalign_path, 
                                                                           strSim_cutoff)
    print('\nHub protein: {}'.format(hub))
    print('    {} of {} hub proteins have been processed. '.format(count_hub, total_hub))
    print('    Calculating structural similarity for protein pairs: ')
    print('         {} + {}'.format(protein_1, protein_2))
    print('         {} of {} interactor pairs have been processed. '.format(count_pairs, total_pairs))
    print('    TM-score for {} and {} is: '.format(protein_1, protein_2))
    print('         {}, {}, {}'.format(TMscore_50, TMscore_30, TMscore_10))
    return TMscore_50, TMscore_30, TMscore_10

# def progress_bar(progress, total, color=colorama.Fore.YELLOW): 
#     percent = 100 * (progress / float(total))
#     bar = '█' * int(percent) + '-' * (100 - int(percent))
#     print(color + f"\r|{bar}| {percent:.2f}%", end="\r")
#     if progress == total: 
#         print(colorama.Fore.GREEN + f"\r|{bar}| {percent:.2f}%", end="\r")

def calculate_TM_score_local_extended(path_1, 
                                      pos_1, 
                                      path_2, 
                                      pos_2, 
                                      protein_1, 
                                      protein_2, 
                                      allProtSeq, 
                                      interactingRes1, 
                                      interactingRes2, 
                                      TempFile_1, 
                                      TempFile_2, 
                                      TMalign_file, 
                                      TMalign_path, 
                                      strSim_cutoff): 
    step = 15
    start_1, end_1 = min(interactingRes1), max(interactingRes1)
    start_2, end_2 = min(interactingRes2), max(interactingRes2)
    res_1 = list(range(start_1, end_1+1))
    res_2 = list(range(start_2, end_2+1))
    TMscore_max_50 = -1
    TMscore_max_30 = -1
    TMscore_max_10 = -1
    i = 0
    print(len(allProtSeq[protein_1]), len(allProtSeq[protein_2]))
    ini_r1_min, ini_r1_max, ini_r2_min, ini_r2_max = min(res_1), len(allProtSeq[protein_1])-max(res_1)+1, min(res_2), len(allProtSeq[protein_2])-max(res_2)+1
    while min(res_1) >= 1 or max(res_1) <= len(allProtSeq[protein_1]) or min(res_2) >= 1 or max(res_2) <= len(allProtSeq[protein_2]): 
        with open(path_1, 'r') as fin, open(TempFile_1, 'w') as fout: 
            for line in fin.readlines(): 
                if line.startswith('ATOM') and line[21] == pos_1 and int(line[22:26]) in res_1: 
                    fout.write(line)
        with open(path_2, 'r') as fin, open(TempFile_2, 'w') as fout: 
            for line in fin.readlines(): 
                if line.startswith('ATOM') and line[21] == pos_2 and int(line[22:26]) in res_2: 
                    fout.write(line)
        cmd = [str(TMalign_path), str(TempFile_1), str(TempFile_2), '-a', 'T']
        summary = subprocess.run(cmd, capture_output = True)
        if i > start_1 - 1: 
            intResPos_1 = [pos for pos in range(start_1, end_1 + 1)]
        else: 
            intResPos_1 = [pos for pos in range(i + 1, i + end_1 - start_1 + 2)]
        if i > start_2 - 1: 
            intResPos_2 = [pos for pos in range(start_2, end_2 + 1)]
        else: 
            intResPos_2 = [pos for pos in range(i + 1, i + end_2 - start_2 + 2)]
        TMscore = retrieve_TM_score(summary.stdout)
        if TMscore != -1: 
            if check_mapping_status_for_interface(summary.stdout, intResPos_1, intResPos_2, 0.5): 
                TMscore_max_50 = max(TMscore, TMscore_max_50)
            if check_mapping_status_for_interface(summary.stdout, intResPos_1, intResPos_2, 0.3): 
                TMscore_max_30 = max(TMscore, TMscore_max_30)
            if check_mapping_status_for_interface(summary.stdout, intResPos_1, intResPos_2, 0.1): 
                TMscore_max_10 = max(TMscore, TMscore_max_10)
        delete_file(TMalign_file)
        delete_file(TempFile_1)
        delete_file(TempFile_2)
        res_1 += list(range(min(res_1)-step, min(res_1)))
        res_1 += list(range(max(res_1)+1, max(res_1)+1+step))
        res_2 += list(range(min(res_2)-step, min(res_2)))
        res_2 += list(range(max(res_2)+1, max(res_2)+1+step))
        # res_1.append(min(res_1)-1) 
        # res_1.append(max(res_1)+1)
        # res_2.append(min(res_2)-1)
        # res_2.append(max(res_2)+1)
        i += step
        print(f'{i}', end='\r')
        # progress_bar(i, max([ini_r1_min, ini_r1_max, ini_r2_min, ini_r2_max]))
    return TMscore_max_50, TMscore_max_30, TMscore_max_10