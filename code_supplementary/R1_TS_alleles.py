import pandas as pd
import numpy as np 
from pathlib import Path
from import_tools import (import_structural_interactome, 
                          import_hub_proteins)
from interface_tools import (get_partner_list,
                             calculate_degree_of_interface_overlap)

def read_GIPS(inPath): 
    result_df = pd.read_table(inPath, header = [0,1], index_col=[0,1])
    result_df = result_df.groupby(result_df.columns.get_level_values(1)).mean()
    result_df = result_df.groupby(result_df.columns.get_level_values(1), axis=1).mean()
    return result_df

def main():
    InPathStrInt = Path('../data/processed/Yeast/model_based/structural_interactome.txt')
    InPathIDMappingFile = Path('../data/external/YEAST_559292_idmapping.dat')
    InPathDegreeStr = Path('../data/processed/aim_1/degree_count_structural_interactome.txt')
    GIPS_matrix = read_GIPS('../data/external/GIPS_matrix/cc_ALL.txt')
    OutPath = Path('/../data/processed/round_1/TS_alleles.txt')
    strInt = import_structural_interactome(InPathStrInt, InPathIDMappingFile)
    hubProteins, hubDict = import_hub_proteins(InPathDegreeStr, InPathIDMappingFile, 2)
    result_all = []
    count = 0
    for hub in hubProteins: 
        partnerList = get_partner_list(strInt, hub)
        strInt_temp = strInt[(strInt['Protein_1'] == hub) | (strInt['Protein_2'] == hub)]
        for partner_1, partner_2 in [(a, b) for idx, a in enumerate(partnerList) for b in partnerList[idx+1:]]: 
            index = '_'.join(sorted([partner_1, partner_2]))
            ORC, ORCR = calculate_degree_of_interface_overlap(strInt_temp, hub, partner_1, partner_2)
            if not ((partner_1 in GIPS_matrix.columns) and (partner_2 in GIPS_matrix.columns)): 
                continue
            GIPS = GIPS_matrix.loc[partner_1, partner_2]
            result_all.append( (hub, partner_1, partner_2, ORC, ORCR, GIPS, hubDict[hub]) )
        count += 1
        print('\r{} of {} hub proteins are processed. '.format(count, len(hubProteins)), end='')
    result_all_df = pd.DataFrame(result_all, columns=['hub', 'partner_1', 'partner_2', 'ORC', 'ORCR', 'GIPS', 'PPI_degree'])
    result_all_df.dropna(inplace=True)
    result_all_df.to_csv(OutPath, sep='\t', index=False, na_rep='NA')
    print(result_all_df[['ORC', 'GIPS']].corr())
    X = result_all_df['ORC'].values.reshape(-1, 1)
    y = result_all_df['GIPS'].values
    X_with_intercept = np.c_[np.ones(X.shape[0]), X]
    coefficients = np.linalg.inv(X_with_intercept.T @ X_with_intercept) @ X_with_intercept.T @ y
    intercept, slope = coefficients
    print(f"Fitted Line: GIPS = {intercept} + {slope} * ORC")

if __name__ == '__main__':
    main()