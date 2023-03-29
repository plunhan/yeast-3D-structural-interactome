import pandas as pd 
import numpy as np
import pickle

def integrate_protein_interaction_data(InputList, IDMappingFile, OutPath): 
    frames = []
    for datasetPath in InputList: 
        dataset = pd.read_csv(datasetPath, sep=' ')[['orf_name_a', 'orf_name_b']]
        dataset['orf_name_a'], dataset['orf_name_b'] = dataset.min(axis=1), dataset.max(axis=1)
        # Create a tag for each protein-protein interaction
        # dataset['index'] = dataset.min(axis=1) + '_' + dataset.max(axis=1) 
        frames.append(dataset)
    output = pd.concat(frames)
    output = output[output['orf_name_a'] != output['orf_name_b']]
    IDMappingDict = {}
    with open(IDMappingFile, 'r') as f:
        lines = f.readlines()
        for line in lines: 
            columns = line.rstrip().split('\t')
            if columns[1] == 'Gene_OrderedLocusName': 
                IDMappingDict[columns[2]] = columns[0]
    output = output[(output['orf_name_a'].isin(IDMappingDict.keys())) & (output['orf_name_b'].isin(IDMappingDict.keys()))]
    output['orf_name_a'] = output['orf_name_a'].map(IDMappingDict)
    output['orf_name_b'] = output['orf_name_b'].map(IDMappingDict)
    output.columns = ['Protein_1', 'Protein_2']
    output.to_csv(OutPath, sep='\t', index=False)

def process_BIOGRID_protein_interactome(InPath, OutPath, species='cerevisiae'): 
    BIOGRID = pd.read_table(InPath)
    if species == 'cerevisiae': 
        species = 'Saccharomyces cerevisiae (S288c)'
        organismID = 559292
    elif species == 'pombe': 
        species = 'Schizosaccharomyces pombe (972h)'
        organismID = 284812
    BIOGRID = BIOGRID[BIOGRID['Experimental System Type'] == 'physical']
    BIOGRID = BIOGRID[(BIOGRID['Organism ID Interactor A'] == organismID) & (BIOGRID['Organism ID Interactor B'] == organismID)]
    BIOGRID = BIOGRID[(BIOGRID['SWISS-PROT Accessions Interactor A'] != '-') & (BIOGRID['SWISS-PROT Accessions Interactor B'] != '-')]
    BIOGRID = BIOGRID[['SWISS-PROT Accessions Interactor A', 'SWISS-PROT Accessions Interactor B']]
    BIOGRID.columns = ['Protein_1', 'Protein_2']
    BIOGRID['Protein_1'], BIOGRID['Protein_2'] = BIOGRID[['Protein_1', 'Protein_2']].min(axis=1), BIOGRID[['Protein_1', 'Protein_2']].max(axis=1)
    BIOGRID = BIOGRID.drop_duplicates()
    BIOGRID.to_csv(OutPath, sep='\t', index=False)

def process_SGA_reliable(InPath, OutPath):
    entryDict = {} # the dict is used to compare different entries for a specific gene pair
    with open(InPath, 'r') as f:
        f.readline()
        while True:
            line = f.readline().strip().split('\t')
            if len(line) < 5:
                break
            query, array, GI, pvalue = line[0].split('_')[0], line[2].split('_')[0], line[5], line[6]
            if float(pvalue) > 0.05:
                continue
            query, array = min(query, array), max(query, array)
            index = '_'.join([query, array])
            if (index in entryDict.keys()) and (float(pvalue) > entryDict[index][3]):
                continue
            entryDict[index] = [query, array, float(GI), float(pvalue)] # take the entry with lowest p-value
    entryDf = pd.DataFrame().from_dict(entryDict, orient='index', columns=['query', 'array', 'GI', 'pvalue']).reset_index()
    entryDf.to_csv(OutPath, sep='\t')

def merge_SGA(InPathList, OutPath):
    DataFrameList = []
    for i in range(len(InPathList)):
        df = pd.read_table(InPathList[i])
        DataFrameList.append(df)
    output = pd.concat(DataFrameList)
    output = output.sort_values(by=['index'])
    output = output.drop('Unnamed: 0', axis=1)
    output.to_csv(OutPath, sep='\t', index=False)

def generate_GI_matrix(InPath, OutPath): 
    GI = pd.read_table(InPath)
    geneList = list(set(GI['query']).union(set(GI['array'])))
    geneList.sort()
    zero_data = np.zeros(shape=(len(geneList), len(geneList)))
    GI_matrix = pd.DataFrame(zero_data, columns=geneList, index=geneList)
    print(GI.columns)
    for index, row in GI.iterrows(): 
        if row['GI'] < -0.08: 
            GI_score = -1
        elif row['GI'] > 0.08: 
            GI_score = 1
        else: 
            GI_score = 0
        GI_matrix.loc[row['query'], row['array']] = GI_score
        GI_matrix.loc[row['array'], row['query']] = GI_score
        if index % 10 == 0: 
            print('{} out of {} queries have been processed. '.format(index, len(GI)))
    GI_matrix.to_csv(OutPath, sep='\t')

def generate_essential_gene_list(OutPath): 
    ## Adapted from YeastMine

    # The following two lines will be needed in every python script:
    from intermine.webservice import Service
    service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")

    # Get a new query on the class (table) you will be querying:
    query = service.new_query("Phenotype")

    # The view specifies the output columns
    query.add_view(
        "genes.secondaryIdentifier", "genes.symbol",
        "genes.sgdAlias", "genes.qualifier", "experimentType", "mutantType",
        "observable", "qualifier", "allele", "alleleDescription",
        "strainBackground", "chemical", "condition", "details", "reporter",
        "publications.pubMedId", "publications.citation"
    )

    # This query's custom sort order is specified below:
    query.add_sort_order("Phenotype.experimentType", "ASC")

    # You can edit the constraint values below
    query.add_constraint("observable", "=", "inviable", code="A")

    # Uncomment and edit the code below to specify your own custom logic:
    # query.set_logic("A")
    essentialGeneList = []
    for row in query.rows():
        essentialGeneList.append(row["genes.secondaryIdentifier"])
    pickle.dump(essentialGeneList, open(OutPath, 'wb'))

# def generate_gene_list(InPath, OutPath): 
#     outputset = set()
#     with open(InPath, 'r') as f: 
#         f.readline()
#         while True: 
#             line = f.readline().strip().split('\t')
#             if len(line) < 5: 
#                 break
#             gene_1 = line[0].split('_')[0]
#             gene_2 = line[2].split('_')[0]
#             outputset.add(gene_1)
#             outputset.add(gene_2)
#     pickle.dump(outputset, open(OutPath, 'wb'))

def generate_GI_full_matrix(InPathList, OutPath): 
    # create an empty matrix to record genetic interaction scores
    entryDict = {}
    for InPath in InPathList: 
        with open(InPath, 'r') as f: 
            flag = 0
            for line in f.readlines():
                if flag == 0: 
                    flag = 1
                    continue
                line = line.strip().split('\t')
                query, array, GI, pvalue = line[0].split('_')[0], line[2].split('_')[0], line[5], line[6]
                query, array = min(query, array), max(query, array)
                index = '_'.join([query, array])
                if (index in entryDict.keys()) and (float(pvalue) > entryDict[index][3]):
                    continue
                entryDict[index] = [query, array, float(GI), float(pvalue)]
    GI = pd.DataFrame().from_dict(entryDict, orient='index', columns=['query', 'array', 'GI', 'pvalue']).reset_index()
    geneList = list(set(GI['query']).union(set(GI['array'])))
    geneList.sort()
    zero_data = np.zeros(shape=(len(geneList), len(geneList)))
    GI_matrix = pd.DataFrame(zero_data, columns=geneList, index=geneList)
    GI_matrix[:] = np.nan
    count = 0
    for index, row in GI.iterrows(): 
        GI_matrix.loc[row['query'], row['array']] = row['GI']
        GI_matrix.loc[row['array'], row['query']] = row['GI']
        count = count + 1
        print('\r{:.2f}%% genetic interaction scores have been processed. '.format(count * 100 / len(GI)), end='')
    GI_matrix.to_csv(OutPath, sep='\t')