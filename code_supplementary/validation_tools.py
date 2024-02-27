import os
import shutil
import pandas as pd
from pathlib import Path
from interactome_tools import read_chain_annotated_interactome, write_chain_annotated_interactome

tp_all = 0
fn_all = 0
fp_all = 0
tn_all = 0
all_analyzed = 0
tpr_ls = []

def read_pdb_IDs(pdbIDFile): 
    with open(pdbIDFile, 'r') as f: 
        lines = f.readlines()
        pdbIDs = set([line.strip() for line in lines])
    return pdbIDs

def filter_pdb_files(InPath, pdbDir, pdbYeastDir, OutPath, speciesTaxID = '559292'): 
    pdb_list = []
    i = 0

    with open(InPath, 'r') as f: 
        lines = f.readlines()
        pdbIDs = [line.strip() for line in lines]
    
    for pdbID in pdbIDs: 
        filePath = pdbDir / ('pdb' + pdbID + '.ent')
        print('\rChecking if {} is a yeast PDB structure. {} of {} structures finished. '.format(pdbID, i, len(pdbIDs)), end='')
        if filePath.is_file() and is_species(filePath, str(speciesTaxID)): 
            pdb_list.append(pdbID)
            pastePath = pdbYeastDir / ('pdb' + pdbID + '.ent')
            shutil.copy(filePath, pastePath)
        i += 1

    with open(OutPath, 'w') as fout: 
        for pdb in pdb_list: 
            fout.write(pdb + '\n')

def is_species(InPath, speciesTaxID = '559292'): 
    with open(InPath, 'r') as f: 
        lines = f.readlines()
        for line in lines: 
            linesplit = line.strip().split()
            record = linesplit[0]
            if record == 'SOURCE': 
                if 'ORGANISM_TAXID:' in linesplit: 
                    taxID_index = linesplit.index('ORGANISM_TAXID:') + 1
                    taxID = linesplit[taxID_index].strip().split(';')[0]
                    if taxID == speciesTaxID: 
                        return True
    return False

def filter_chain_annotations_by_template(InPath, pdbIDFile, OutPath): 
    pdbIDs = read_pdb_IDs(pdbIDFile)

    with open(InPath, 'r') as f, open(OutPath, 'w') as fout: 
        headerLine = f.readline()
        fout.write(headerLine)
        header = headerLine.strip().split()
        colChain = header.index('Subject')
        colQlen = header.index('Qlen')
        colSlen = header.index('Slen')
        colIdentities = header.index('Identities')
        lines = f.readlines()
        for line in lines: 
            items = line.strip().split()
            chainID = items[colChain]
            pdbID = chainID.split('_')[0]
            if pdbID in pdbIDs: 
                Qlen = items[colQlen]
                Slen = items[colSlen]
                Identities = items[colIdentities]
                if is_crystal(Qlen, Slen, Identities): 
                    fout.write(line)

def is_crystal(Qlen, Slen, Identities):
    if Qlen == Identities and Slen == Identities: 
        return True
    return False

def filter_chain_annotated_interactome_by_template(InPath, pdbIDFile, OutPath): 
    pdbIDs = read_pdb_IDs(pdbIDFile)

    chainAnnotInt = read_chain_annotated_interactome(InPath)
    chainAnnotInt['Mapping_chains'] = chainAnnotInt.apply(lambda x:
                                    filter_yeast_template_chain_pairs(x['Mapping_chains'], 
                                                                      pdbIDs), 
                                    axis=1)
    chainAnnotInt = chainAnnotInt[chainAnnotInt['Mapping_chains'].apply(len) > 0].reset_index(drop=True)
    write_chain_annotated_interactome(chainAnnotInt, OutPath)

def filter_yeast_template_chain_pairs(mapping_chains, pdbIDs): 
    is_yeast_template = []
    for i, chain_tuple in enumerate(mapping_chains): 
        if chain_tuple[0].split('_')[0] in pdbIDs: 
            is_yeast_template.append(i)
    return [mapping_chains[i] for i in is_yeast_template]

def filter_chain_annotations_by_cocrystal_protein (inPath, proteins, outPath, identical = True):
    """Filter protein chain annotations by proteins.

    Args:
        inPath (Path): path to tab-deleimited file containing protein-chain alignments.
        proteins (list): proteins to be selected.
        outPath (Path): file path to save filtered alignments to.

    """
    with open(inPath, "r", errors='ignore') as f, open(outPath, "w") as fout:
        headers = f.readline()
        fout.write(headers)
        headerSplit = headers.strip().split('\t')
        numCol = len(headerSplit)
        colQuery = headerSplit.index('Query')
        colChain = headerSplit.index('Subject')
        colQlen = headerSplit.index('Qlen')
        colSlen = headerSplit.index('Slen')
        colIdentities = headerSplit.index('Identities')
    with open(inPath, "r", errors='ignore') as f, open(outPath, "a") as fout:
        next(f)
        for line in f:
            items = line.strip().split('\t')
            if len(items) == numCol:
                Qlen = items[colQlen]
                Slen = items[colSlen]
                Identities = items[colIdentities]
                if items[colQuery] in proteins and is_crystal(Qlen, Slen, Identities):
                    fout.write(line)

def jaccard_similarity(list_1, list_2): 
    return len(set(list_1).intersection(set(list_2))) / len(set(list_1).union(set(list_2)))

def generate_ppi_index(df): 
    df['p1'], df['p2'] = df[['Protein_1', 'Protein_2']].min(axis=1), df[['Protein_1', 'Protein_2']].max(axis=1)
    df['index'] = df[['p1', 'p2']].agg('_'.join, axis=1)
    return df

def same_template(template_1, template_2): 
    if sorted(template_1.split('-')) == sorted(template_2.split('-')): 
        return True
    return False

def interface_str_to_list(interface_str): 
    interface_1, interface_2 = interface_str.split('+')
    interface_1, interface_2 = interface_1.split(','), interface_2.split(',')
    return interface_1, interface_2

def calculate_jaccard_similarity_entry(interfaces_1, interfaces_2): 
    i11, i12 = interface_str_to_list(interfaces_1)
    i21, i22 = interface_str_to_list(interfaces_2)
    return jaccard_similarity(i11, i21), jaccard_similarity(i12, i22)

def chain_pairs(template): 
    structure, chain1, chain2 = template.split('-')
    c1 = '_'.join([structure, chain1])
    c2 = '_'.join([structure, chain2])
    return ['+'.join([c1, c2]), '+'.join([c2, c1])]

def load_structural_interactome(path, template=False): 
    df = pd.read_table(path)
    df = generate_ppi_index(df)
    if template: 
        df['Chain_pairs'] = df['Chain_pairs'].str.split('|')
    else: 
        df['Template_file_ID'] = df.apply(lambda x: x['Template_file_ID'].replace('!', ''), axis=1)
    return df

def calculate_jaccard_similarity(InPathGoldStandard, InPathModelled, InPathTemplate, OutPath): 
    gs = load_structural_interactome(InPathGoldStandard)
    md = load_structural_interactome(InPathModelled)
    tp = load_structural_interactome(InPathTemplate, template=True)
    md = md[md['index'].isin(set(gs['index']))]
    gs = gs[gs['index'].isin(set(md['index']))]
    md = md.sort_values(by='index')
    gs = gs.sort_values(by='index')
    result = []
    for i in range(len(gs)): 
        if same_template(gs.iloc[i]['Template_file_ID'], md.iloc[i]['Template_file_ID']): 
            continue
        gsChains = tp.iloc[i]['Chain_pairs']
        mdChains = chain_pairs(md.iloc[i]['Template_file_ID'])
        if any([i in gsChains for i in mdChains]): 
            continue
        j1, j2 = calculate_jaccard_similarity_entry(gs.iloc[i]['Interfaces'], md.iloc[i]['Interfaces'])
        result.append( (gs.iloc[i]['Protein_1'], gs.iloc[i]['Protein_2'], j1, j2) )
    resultdf = pd.DataFrame(result, columns = ['Protein_1', 'Protein_2', 'Jaccard_1', 'Jaccard_2'])
    resultdf.to_csv(OutPath, sep='\t', index=False)

def calculate_statistic_numbers_PPI(gs_interfaces, 
                                    md_interfaces, 
                                    universal_1, 
                                    universal_2, 
                                    index, 
                                    protein_1, 
                                    protein_2, 
                                    template_file_id, 
                                    alignment_file_id, 
                                    tempAlignDict): 
    global tp_all, fn_all, tn_all, fp_all, tpr_ls
    gs_interface_1, gs_interface_2 = interface_str_to_list(gs_interfaces)
    md_interface_1, md_interface_2 = interface_str_to_list(md_interfaces)
    tp1, fn1 = calculate_statistic_numbers_protein(gs_interface_1, md_interface_1, universal_1)
    tp2, fn2 = calculate_statistic_numbers_protein(gs_interface_2, md_interface_2, universal_2)
    complex0, templateID = alignment_file_id.split('_')
    templateID = templateID.replace('!', '')
    p1, p2 = complex0.split('=')
    ch1, ch2 = templateID.split('-')[1], templateID.split('-')[2]
    Evalue1 = tempAlignDict[(p1, '_'.join([template_file_id, ch1]))]
    Evalue2 = tempAlignDict[(p2, '_'.join([template_file_id, ch2]))]
    # print([protein_1, p1, protein_2, p2])
    if protein_1 != p1: 
        Evalue1, Evalue2 = Evalue2, Evalue1
    tpr_ls.append(( protein_1, protein_2, alignment_file_id, Evalue1, Evalue2, index, (tp1+tp2) / (tp1+tp2+fn1+fn2) ) )

def calculate_statistic_numbers_protein(gs_interface, md_interface, universal): 
    global tp_all, fn_all, tn_all, fp_all
    tp = len([r for r in md_interface if r in gs_interface])
    tp_all += tp
    fn = len([r for r in gs_interface if r not in md_interface])
    fn_all += fn
    fp = len([r for r in md_interface if r not in gs_interface])
    fp_all += fp
    tn = universal - tp - fn - fp
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

def reset_statistical_numbers(): 
    global tp_all, tn_all, fp_all, fn_all, all_analyzed
    tp_all = 0
    tn_all = 0
    fp_all = 0
    fn_all = 0
    all_analyzed = 0

def generate_template_alignment_dict(tempAlign): 
    result_dict = tempAlign.set_index(['Query', 'Subject'])['Expect'].to_dict()
    return result_dict

def calculate_positive_rates(InPathGoldStandard, 
                             InPathModelled, 
                             InPathTemplate, 
                             InPathChainMap, 
                             ppiTempAlignmentFile, 
                             OutPath): 
    reset_statistical_numbers()
    global tp_all, tn_all, fp_all, fn_all, all_analyzed, tpr_ls
    gs = load_structural_interactome(InPathGoldStandard)
    md = load_structural_interactome(InPathModelled)
    tp = load_structural_interactome(InPathTemplate, template=True)
    tempAlign = pd.read_table(ppiTempAlignmentFile)
    tempAlignDict = generate_template_alignment_dict(tempAlign)
    proteinLength = read_protein_length_dict(InPathChainMap)
    md = md[md['index'].isin(set(gs['index']))]
    gs = gs[gs['index'].isin(set(md['index']))]
    tp = tp[tp['index'].isin(set(gs['index']))]
    md = md.sort_values(by='index')
    gs = gs.sort_values(by='index')
    tp = tp.sort_values(by='index')
    result = []
    for i in range(len(gs)): 
        index = '_'.join(sorted([gs['Protein_1'].iloc[i], gs['Protein_2'].iloc[i]]))
        if same_template(gs.iloc[i]['Template_file_ID'], md.iloc[i]['Template_file_ID']): 
            continue
        gsChains = tp.iloc[i]['Chain_pairs']
        mdChains = chain_pairs(md.iloc[i]['Template_file_ID'])
        if any([i in gsChains for i in mdChains]): 
            continue
        calculate_statistic_numbers_PPI(gs.iloc[i]['Interfaces'], 
                                        md.iloc[i]['Interfaces'], 
                                        proteinLength[gs.iloc[i]['Protein_1']], 
                                        proteinLength[gs.iloc[i]['Protein_2']], 
                                        index, 
                                        gs.iloc[i]['Protein_1'], 
                                        gs.iloc[i]['Protein_2'], 
                                        md.iloc[i]['Template_file_ID'], 
                                        md.iloc[i]['Alignment_file_ID'], 
                                        tempAlignDict)
        all_analyzed += 1
    tpr = tp_all / (tp_all + fn_all)
    fpr = fp_all / (fp_all + tn_all)
    # tpr_ls.append(( protein_1, protein_2, alignment_file_id, Evalue1, Evalue2, index, (tp1+tp2) / (tp1+tp2+fn1+fn2) ) )
    tpr_df = pd.DataFrame(tpr_ls, columns=['Protein_1', 'Protein_2', 'Alignment_file_ID', 'Evalue1', 'Evalue2', 'PPI', 'TPR'])
    if not OutPath.is_file(): 
        tpr_df.to_csv(OutPath, sep='\t', index=False)
    print("In total, there are {} PPIs. ".format(all_analyzed))
    print("True positive prediction is {}. False negative prediction is {}. True positive rate is {}".format(tp_all, fn_all, tpr))
    print("False positive prediction is {}. True negative prediction is {}. False positive rate is: {}".format(fp_all, tn_all, fpr))

def check_PPI_cocrystal(protein_1, protein_2, chain_1, chain_2, cocrystalDict): 
    pair_1 = '+'.join([protein_1, chain_1])
    pair_2 = '+'.join([protein_2, chain_2])
    if pair_1 in cocrystalDict.keys() and pair_2 in cocrystalDict.keys(): 
        return True
    return False

def check_number_of_cocrystal_PPIs(InpathStrInt, InPathChainMap): 
    count = 0
    strInt = load_structural_interactome(InpathStrInt)
    cocrystalDict = {}
    with open(InPathChainMap, 'r') as f: 
        headers = f.readline()
        headerSplit = headers.strip().split('\t')
        print(headerSplit)
        numCol = len(headerSplit)
        colQuery = headerSplit.index('Query')
        colChain = headerSplit.index('Subject')
        colQlen = headerSplit.index('Qlen')
        colSlen = headerSplit.index('Slen')
        colIdentities = headerSplit.index('Identities')
        for line in f.readlines(): 
            items = line.strip().split('\t')
            if len(items) == numCol:
                protein = items[colQuery]
                Qlen = items[colQlen]
                chain = items[colChain]
                Slen = items[colSlen]
                Identities = items[colIdentities]
                if is_crystal(Qlen, Slen, Identities):
                    cocrystalDict['+'.join([protein, chain])] = 1
    for index, row in strInt.iterrows(): 
        proteins, template = row['Alignment_file_ID'].split('_')
        protein_1, protein_2 = proteins.split('=')
        template_id, chainID_1, chainID_2 = template.split('-')
        chain_1, chain_2 = '_'.join([template_id, chainID_1]), '_'.join([template_id, chainID_2])
        if check_PPI_cocrystal(protein_1, protein_2, chain_1, chain_2, cocrystalDict): 
            count += 1
    print('{} of {} PPIs are co-crystal. '.format(count, len(strInt)))