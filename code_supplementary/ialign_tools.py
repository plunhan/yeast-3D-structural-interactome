import pandas as pd
import subprocess
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

def model_name(protein_1, protein_2, models): 
    model = '='.join([protein_1, protein_2])
    if model not in models: 
        model = '='.join([protein_2, protein_1])
    return 'pdb' + model + '.ent'

def calculate_interfacial_similarity(ialign, 
                                     hub, 
                                     protein_1, 
                                     protein_2, 
                                     models, 
                                     ialignTempFile, 
                                     pdbDir, 
                                     distance=5, 
                                     method='tm'): 
    model_1 = model_name(hub, protein_1, models)
    model_2 = model_name(hub, protein_2, models)
    cmd = ['perl', str(ialign), '-a', '2', '-e', method, '-o', str(ialignTempFile), 
           '-dc', str(distance), str(pdbDir / model_1), str(pdbDir / model_2)]
    print(' '.join(cmd))
    print('Running ialign...')
    subprocess.run(cmd, shell=True)
    print('Finished. ')
    intStrSim = retrieve_interfacial_similarity(ialignTempFile, method)
    return intStrSim

def retrieve_interfacial_similarity(ialignFile, method = 'tm'): 
    if method == 'tm': 
        key = 'TM-score'
    elif method == 'is': 
        key = 'IS-score'
    else: 
        print('"method" can only be set to either "tm" or "is". ')
        return False
    with open(ialignFile, 'r') as f: 
        for line in f.readlines(): 
            if key in line: 
                score = float(line.split(' ')[2].split(',')[0])
                return score
    print('Cannot find the score under the method. ')
    return False

def test_interfacial_similarity(strIntPath, 
                                degreeStrIntPath, 
                                IDMappingFile, 
                                GImatrixPath, 
                                ialignPath, 
                                ialignTempFile, 
                                pdbDir, 
                                outPath, 
                                hub_cutoff=2,
                                distance=5, 
                                method='tm'): 
    strInt = pd.read_csv(strIntPath, sep='\t')
    models = set(list(strInt['Complex_ID']))
    degreeStrInt = pd.read_csv(degreeStrIntPath, sep='\t')
    GI_profile = import_GI_profile(GImatrixPath)
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    hubProteins = sorted(list(set(degreeStrInt[degreeStrInt['Degree'] >= hub_cutoff]['Protein'])))
    i = 1
    l = len(hubProteins)
    for hub in hubProteins: 
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        partnerList = get_partner_list(strInt_temp, hub)
        print(f'Finished finding partners for {i} in {l} hub proteins. ')
        partnerPairList = [(partner_1, partner_2) for idx, partner_1 in enumerate(partnerList) for partner_2 in partnerList[idx+1:]]
        for partnerPair in partnerPairList: 
            protein_1, protein_2 = partnerPair
            if protein_1 == protein_2: 
                continue
            print(f'Calculating interfacial structural similarity for {protein_1} and {protein_2}. ')
            intStrSim = calculate_interfacial_similarity(ialignPath, 
                                                         hub, 
                                                         protein_1, 
                                                         protein_2, 
                                                         models, 
                                                         ialignTempFile, 
                                                         pdbDir, 
                                                         distance=distance, 
                                                         method=method)
            print(f'Calculating interfacial overlap for {protein_1} and {protein_2}. ')
            if not intStrSim: 
                continue
            ORC = calculate_degree_of_interface_overlap(strInt_temp, hub, protein_1, protein_2)
            protein_1, protein_2 = IDMappingDict[protein_1], IDMappingDict[protein_2]
            print(f'Calculating genetic interaction profile similarity for {protein_1} and {protein_2}. ')
            GIPS = calculate_GIPS(GI_profile, protein_1, protein_2)
            if not GIPS or np.isnan(GIPS): 
                continue
            resultList.append( (IDMappingDict[hub], protein_1, protein_2, ORC, GIPS, intStrSim) )
        i += 1
    resultDataFrame = pd.DataFrame(resultList, columns = ['hub', 'Protein_1', 'Protein_2', 'ORC', 'GIPS', 'TM-score'])
    resultDataFrame.to_csv(outPath, sep = '\t', na_rep = 'NA')