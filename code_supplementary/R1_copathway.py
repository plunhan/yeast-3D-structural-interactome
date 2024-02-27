import os
from pathlib import Path
import pandas as pd
from import_tools import map_protein_id_to_locus_id

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

    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)

    pathbankFile = extDir / 'pathbank_all_proteins.csv'
    pathbank = pd.read_csv(pathbankFile)
    pb_yeast = pathbank[pathbank['Species'] == 'Saccharomyces cerevisiae']
    pb_dict = dict()

    for index, row in pb_yeast.iterrows(): 
        protein = row['Uniprot ID']
        pathway = row['PathBank ID']
        pathwayName = row['Pathway Name']
        if protein not in IDMappingDict.keys(): 
            continue
        if protein not in pb_dict.keys(): 
            pb_dict[IDMappingDict[protein]] = [(pathway, pathwayName)]
        else: 
            pb_dict[IDMappingDict[protein]].append( (pathway, pathwayName) )

    def get_pathway(hub, partner_1, partner_2, pb_dict): 
        if any(protein not in pb_dict.keys() for protein in [hub, partner_1, partner_2]): 
            return 'Not Applicable'
        s1 = set(pb_dict[hub])
        s2 = set(pb_dict[partner_1])
        s3 = set(pb_dict[partner_2])
        result_ls = list(s1.intersection(s2, s3))
        if len(result_ls): 
            return ','.join(sorted([result[1] for result in result_ls]))
        return 'Not Applicable'

    corResultFile = aim2Dir / 'overlapping_residue_count_vs_genetic_GIPS.txt'
    corResult = pd.read_table(corResultFile)
    corResult = corResult.dropna()
    corResult['pathways'] = corResult.apply(lambda row: get_pathway(row['hub'], row['partner_1'], row['partner_2'], pb_dict), axis=1)
    print(corResult[corResult['pathways']!='Not Applicable'])

    output = rd1Dir / 'result_copathway.txt'
    corResult[corResult['pathways']!='Not Applicable'].to_csv(output, sep='\t', index=False)

if __name__ == '__main__':
    main()