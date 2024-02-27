import pandas as pd 
from pathlib import Path 
import networkx as nx

def map_protein_id_to_locus_id(IDMappingFile): 
    IDMappingDict = {}
    with open(IDMappingFile, 'r') as f: 
        lines = f.readlines()
        for line in lines: 
            elements = line.rstrip().split('\t')
            if elements[1] == 'Gene_OrderedLocusName': 
                IDMappingDict[elements[0]] = elements[2]
    return IDMappingDict

def main(): 
    dataDir = Path('../data/processed/round_1/bin_tmscore')
    dataframe = dataDir / 'bin_tmscore.txt'
    resultPath = Path('../data/processed/aim_2/overlapping_residue_count_vs_genetic_GIPS.txt')
    result = pd.read_table(resultPath)
    IDMappingDict = map_protein_id_to_locus_id('../data/external/YEAST_559292_idmapping.dat')
    G = nx.Graph()
    df = pd.read_table(dataframe)
    df.columns = ['p1', 'p2', 'tm']
    for row in df.itertuples(index=False): 
        if row.tm > 0.5: 
            G.add_edge(row.p1, row.p2, tmscore=row.tm)
        else: 
            G.add_node(row.p1)
            G.add_node(row.p2)

    group_labels = {IDMappingDict[hub]: idx for idx, hubs in enumerate(nx.connected_components(G)) for hub in hubs}
    result['group'] = result['hub'].map(group_labels)
    result = result.groupby('group').agg({'ORC': 'mean', 'GIPS': 'mean'}).reset_index()
    result.to_csv(dataDir / 'binned_ORC_GIPS.txt', sep='\t', index=False)

if __name__ == '__main__':
    main()