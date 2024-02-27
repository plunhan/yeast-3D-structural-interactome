import os
from pathlib import Path
from import_tools import map_protein_id_to_locus_id
from itertools import combinations
from scipy.stats import pearsonr
import subprocess
import re
import pandas as pd
import multiprocessing as mp

def delete_file(file): 
    if file.is_file(): 
        subprocess.call(['rm', str(file)])

def retrieve_model_name(hub, partner, ppiModelDir): 
    file_1 = 'pdb' + hub + '=' + partner + '.ent'
    file_2 = 'pdb' + partner + '=' + hub + '.ent'
    path_1 = ppiModelDir / file_1
    path_2 = ppiModelDir / file_2
    if path_1.is_file(): 
        path = path_1
        hub_pos = 'A'
    elif path_2.is_file(): 
        path = path_2
        hub_pos = 'B'
    else: 
        return '-1', '-1'
    return (path, hub_pos)

def create_model_dict(cor, ppiModelDir): 
    cor['model_name'] = cor.apply(lambda row: retrieve_model_name(row['hub'], row['partner_1'], ppiModelDir), axis=1)
    model_dict = dict(zip(cor['hub'], cor['model_name']))
    return model_dict

def retrieve_TM_score(line): 
    try: 
        TMscore = float(re.findall(r"\b\d+\.\d+\s\(", str(line))[2].split(' ')[0])
        return TMscore
    except: 
        return -1

def calculate_hub_structural_similarity(hub1, 
                                        hub2, 
                                        model_dict, 
                                        TempFile_1, 
                                        TempFile_2, 
                                        TMalign_file, 
                                        TMalign_path): 
    path_1, pos_1 = model_dict[hub1]
    path_2, pos_2 = model_dict[hub2]
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

def bin_tmscore_worker(args):
    hub1, hub2, model_dict, TempFile_1, TempFile_2, TMalign_file, TMalign_path = args
    TempFile_1 = Path(str(TempFile_1).split('.txt')[0] + '_' + hub1 + '_' + hub2 + '1.txt')
    TempFile_2 = Path(str(TempFile_1).split('.txt')[0] + '_' + hub1 + '_' + hub2 + '2.txt')
    TMscore = calculate_hub_structural_similarity(hub1, hub2, model_dict, TempFile_1, TempFile_2, TMalign_file, TMalign_path)
    return (hub1, hub2, TMscore)

def bin_tmscore(corResult, IDMappingFile, ppiModelDir, TempFile_1, TempFile_2, TMalign_file, TMalign_path, binTMscoreResult):
    ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
    pool = mp.Pool(processes=ncpus)

    cor = pd.read_table(corResult)
    cor.dropna(inplace=True)
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    IDMappingDict = {value: key for key, value in IDMappingDict.items()}
    cor.replace(IDMappingDict, inplace=True)
    cor.drop_duplicates(subset='hub', keep='first', inplace=True)
    hub_ls = sorted(list(set(cor['hub'])))
    model_dict = create_model_dict(cor, ppiModelDir)
    friend_circles = dict()
    all_list = list(combinations(sorted(hub_ls), 2))
    all_list = all_list[:1000]

    worker_args = [(hub1, hub2, model_dict, TempFile_1, TempFile_2, TMalign_file, TMalign_path) for hub1, hub2 in all_list]

    result = pool.map(bin_tmscore_worker, worker_args)
    result_df = pd.DataFrame(result, columns = ['protein_1', 'protein_2', 'TM-score'])
    result_df.to_csv(binTMscoreResult, sep='\t', index=False)

    pool.close()
    pool.join()

def main(): 

    # parent directory of all data files
    dataDir = Path('/home/han1204/scratch')

    # directory of TM-align files
    TMalignDir = dataDir / 'TMalign'

    # input files
    corResult = dataDir / 'overlapping_residue_count_vs_genetic_GIPS.txt'
    IDMappingFile = dataDir / 'YEAST_559292_idmapping.dat'
    TempFile_1 = dataDir / 'tempfile1.txt'
    TempFile_2 = datair / 'tempfile2.txt'
    TMalign_file = dataDir / 'TMscore.txt'

    # output files
    binTMscoreResult = dataDir / 'bin_tmscore.txt'

    if not binTMscoreResult.is_file(): 
        bin_tmscore(corResult, 
                    IDMappingFile, 
                    ppiModelDir, 
                    TempFile_1, 
                    TempFile_2, 
                    TMalign_file, 
                    TMalignDir, 
                    binTMscoreResult)

if __name__ == '__main__':
    main()