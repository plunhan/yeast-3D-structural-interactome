import pandas as pd
from pathlib import Path

tp_all = 0
tn_all = 0
fp_all = 0
fn_all = 0
all_analyzed = 0
tpr_ls = []

def interface_str_to_list(interface_str): 
    interface_1, interface_2 = interface_str.split('+')
    interface_1, interface_2 = interface_1.split(','), interface_2.split(',')
    return interface_1, interface_2

def calculate_statistic_numbers_protein(gs_interface, md_interface, protein_len): 
    global tp_all, fn_all, tn_all, fp_all
    tp = len([r for r in md_interface if r in gs_interface])
    tp_all += tp
    fn = len([r for r in gs_interface if r not in md_interface])
    fn_all += fn
    fp = len([r for r in md_interface if r not in gs_interface])
    fp_all += fp
    tn = protein_len - tp - fn - fp
    tn_all += tn
    return tp, fn 

def read_protein_length_dict(InPath): 
    proteinLengthDict = {}
    with open(InPath, 'r') as f: 
        headerLine = f.readline()
        headerSplit = headerLine.strip().split('\t')
        colQuery = headerSplit.index('Query')
        colQlen = headerSplit.index('Qlen')
        for line in f.readlines(): 
            linesplit = line.strip().split('\t')
            Query = linesplit[colQuery]
            Qlen = linesplit[colQlen]
            proteinLengthDict[Query] = int(Qlen)
    return proteinLengthDict

def main(): 
    gs = pd.read_table('../data/processed/Yeast/validation/model_based/structural_interactome.txt')
    ap = pd.read_table('../data/processed/round_1/structural_interactome_alphafold.txt', index_col=0)
    OutPath = Path('~/Desktop/RealScience/data/processed/round_1/distribution_positive_rate_alphafold_corrected.txt')
    gs['index'] = gs.apply(lambda row: '_'.join(sorted([row['Protein_1'], row['Protein_2']])), axis=1)
    ap['index'] = ap.apply(lambda row: '_'.join(sorted([row['Protein_1'], row['Protein_2']])), axis=1)
    intersect = set(ap['index']).intersection(set(gs['index']))
    gs = gs[gs['index'].isin(intersect)]
    ap = ap[ap['index'].isin(intersect)]
    gs = gs.sort_values(by='index')
    ap = ap.sort_values(by='index')
    print(gs[['Protein_1', 'Protein_2', 'Interfaces']])
    print(ap[['Protein_1', 'Protein_2', 'Interfaces']])
    InPathChainMap = Path('/Users/plunhan/Desktop/RealScience/data/processed/Yeast/model_based/struc_interactome_chain_map.txt')
    proteinLength = read_protein_length_dict(InPathChainMap)
    for i in range(len(gs)): 
        index = gs['index'].iloc[i]
        p_1, p_2 = index.split('_')
        gs_interface_1, gs_interface_2 = interface_str_to_list(gs['Interfaces'].iloc[i])
        if gs['Protein_1'].iloc[i] == p_2: 
            gs_interface_2, gs_interface_1 = gs_interface_1, gs_interface_2
        ap_interface_1, ap_interface_2 = interface_str_to_list(ap['Interfaces'].iloc[i])
        if ap['Protein_1'].iloc[i] == p_2: 
            ap_interface_2, ap_interface_1 = ap_interface_1, ap_interface_2
        tp1, fn1 = calculate_statistic_numbers_protein(gs_interface_1, ap_interface_1, proteinLength[gs['Protein_1'].iloc[i]])
        tp2, fn2 = calculate_statistic_numbers_protein(gs_interface_2, ap_interface_2, proteinLength[gs['Protein_2'].iloc[i]])
        tpr_ls.append(( index, (tp1+tp2) / (tp1+tp2+fn1+fn2) ))
    tpr = tp_all / (tp_all + fn_all)
    fpr = fp_all / (fp_all + tn_all)
    tpr_df = pd.DataFrame(tpr_ls, columns=['PPI', 'TPR'])
    if not OutPath.is_file(): 
        tpr_df.to_csv(OutPath, sep='\t', index=False)
    print("True positive prediction is {}. False negative prediction is {}. True positive rate is {}".format(tp_all, fn_all, tpr))
    print("False positive prediction is {}. True negative prediction is {}. False positive rate is: {}".format(fp_all, tn_all, fpr))

if __name__ == '__main__':
    main()