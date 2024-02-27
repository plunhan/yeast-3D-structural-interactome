import pandas as pd
from pathlib import Path
import os
from scipy.stats import pearsonr
from import_tools import map_protein_id_to_locus_id

def main():
    # input files
    corResult = '../data/processed/aim_2/overlapping_residue_count_vs_genetic_GIPS.txt'
    complexomeFile = '../data/external/559292.tsv'
    IDMappingFile = '../data/external/YEAST_559292_idmapping.dat'
    # output files
    binnedComplex = '../data/processed/round_1/bin_complex.txt'

    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    complexome = pd.read_csv(complexomeFile, sep = '\t')
    complexome = complexome[['#Complex ac', 'Identifiers (and stoichiometry) of molecules in complex']]
    complexome.columns = ['Complex_name', 'Subunits']

    ternary = dict()
    list_of_lists = []

    for index, row in complexome.iterrows(): 
        complexName = row['Complex_name']
        proteins_stoc = row['Subunits'].split('|')
        proteins = [entry.split('(')[0] for entry in proteins_stoc]
        protein_ls = [IDMappingDict[protein] for protein in proteins if protein in IDMappingDict.keys()]
        protein_ls.sort()
        ternary[complexName] = protein_ls
        list_of_lists.append(protein_ls)

    def check_strings_occurrence(tuple_value, ternary):
        for complexName, protein_ls in ternary.items(): 
            if all(string in protein_ls for string in tuple_value): 
                return complexName
        return '_'.join(tuple_value)

    def create_tuple(row): 
        return (row['hub'], row['partner_1'], row['partner_2'])

    overlap = pd.read_csv(corResult, sep='\t')
    overlap.dropna(inplace=True)
    overlap['index'] = overlap.apply(create_tuple, axis=1)
    overlap['complex'] = overlap['index'].apply(lambda x: check_strings_occurrence(x, ternary))
    overlap_ternary = overlap.groupby('complex').mean().reset_index()
    overlap_ternary.to_csv(binnedComplex, sep='\t', index=False)
    print(pearsonr(overlap_ternary['ORC'], overlap_ternary['GIPS']))

if __name__ == '__main__':
    main()