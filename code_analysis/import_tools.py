# Import processed data tools 
import pandas as pd 
import pickle

def map_protein_id_to_locus_id(IDMappingFile): 
    '''

    This function returns a dict that map Uniprot ID to gene ordered locus name. 

    Input: 

        IDMappingFile (Path): path to ID mapping file

    Return: 

        dict: key is Uniprot ID and value is gene ordered locus name
    
    '''
    IDMappingDict = {}
    with open(IDMappingFile, 'r') as f: 
        lines = f.readlines()
        for line in lines: 
            elements = line.rstrip().split('\t')
            if elements[1] == 'Gene_OrderedLocusName': 
                IDMappingDict[elements[0]] = elements[2]
    return IDMappingDict

def map_protein_id_to_gene_name(IDMappingFile): 
    '''

    This function returns a dict that map Uniprot ID to gene name. 

    Input: 

        IDMappingFile (Path): path to ID mapping file

    Return: 

        dict: key is Uniprot ID and value is gene ordered locus name
    
    '''
    IDMappingDict = {}
    with open(IDMappingFile, 'r') as f: 
        lines = f.readlines()
        for line in lines: 
            elements = line.rstrip().split('\t')
            if elements[1] == 'Gene_Name': 
                IDMappingDict[elements[0]] = elements[2]
    return IDMappingDict

def map_protein_id_to_ORF_name_pombe(IDMappingFile): 
    IDMappingDict = {}
    with open(IDMappingFile, 'r') as f: 
        lines = f.readlines()
        for line in lines: 
            elements = line.rstrip().split('\t')
            if elements[1] == 'Gene_ORFName': 
                IDMappingDict[elements[0]] = elements[2].upper()
    return IDMappingDict

def import_structural_interactome(InPathStrInt, IDMappingFile, index=False): 
    '''

    This function returns a dataframe that contains interactors' names and interfacial residues.  

    Input: 
        
        InPathStrInt (Path): path to structural interactome
        IDMappingFile (Path): path to ID mapping file

    Return: 

        DataFrame: structural interactome with protein names mapped to gene ordered locus name and without duplicates
    
    '''
    strInt = pd.read_table(InPathStrInt) 
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    strInt['Protein_1'] = strInt['Protein_1'].map(IDMappingDict)
    strInt['Protein_2'] = strInt['Protein_2'].map(IDMappingDict)
    strInt = strInt[strInt['Protein_1'] != strInt['Protein_2']]
    strInt = strInt.drop_duplicates()
    if index: 
        strInt['index'] = strInt[['Protein_1', 'Protein_2']].agg('_'.join, axis=1)
        strInt = strInt.set_index('index')
    return strInt

def import_structural_interactome_with_interacting_residues(InPathStrInt, IDMappingFile): 
    '''

    This function returns a dataframe that contains interactors' names and interacting residues on interface. 

    Input: 
        
        InPathStrInt (Path): path to structural interactome with information of interacting residues
        IDMappingFile (Path): path to ID mapping file

    Return: 

        DataFrame: structural interactome with protein names mapped to gene ordered locus name and without duplicates
    
    '''
    strInt = pickle.load(open(InPathStrInt, 'rb'))
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    strInt['Protein_1'] = strInt['Protein_1'].map(IDMappingDict)
    strInt['Protein_2'] = strInt['Protein_2'].map(IDMappingDict)
    strInt = strInt[strInt['Protein_1'] != strInt['Protein_2']]
    return strInt

def import_hub_proteins(InPathDegree, IDMappingFile, hubCutOff=2): 
    '''

    This function returns a list of hub proteins defined by the cut-off in Aim_1.py. 

    Input: 
        
        InPathDegree (Path): path to the file contains information of degree of PPI
        IDMappingFile (Path): path to ID mapping file

    Return: 

        list: hub proteins defined by the cut-off in Aim_1.py 
    
    '''
    degreeStrInt = pd.read_table(InPathDegree)
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    hubProteins = list(set(degreeStrInt[degreeStrInt['Degree'] >= hubCutOff]['Protein']))
    hubProteins = [IDMappingDict[x] for x in hubProteins if x in IDMappingDict.keys()]
    hubDict = dict(zip(degreeStrInt['Protein'].map(IDMappingDict), degreeStrInt['Degree']))
    return hubProteins, hubDict

def import_human_hub_proteins(InPathDegree, IDMappingFile, hubCutOff=2): 
    '''

    This function returns a list of hub proteins defined by the cut-off in Aim_1.py. 

    Input: 
        
        InPathDegree (Path): path to the file contains information of degree of PPI
        IDMappingFile (Path): path to ID mapping file

    Return: 

        list: hub proteins defined by the cut-off in Aim_1.py 
    
    '''
    degreeStrInt = pd.read_table(InPathDegree)
    IDMappingDict = map_protein_id_to_gene_name(IDMappingFile)
    hubProteins = list(set(degreeStrInt[degreeStrInt['Degree'] >= hubCutOff]['Protein']))
    hubProteins = [IDMappingDict[x] for x in hubProteins if x in IDMappingDict.keys()]
    hubDict = dict(zip(degreeStrInt['Protein'].map(IDMappingDict), degreeStrInt['Degree']))
    return hubProteins, hubDict

def import_sequence_similarity(InPathSeqSim): 
    '''

    This function returns two dicts that record identities and positives for interactor pairs, respectively. 

    Input: 
        
        InPathDegree (Path): path to the file contains information of degree of PPI
        IDMappingFile (Path): path to ID mapping file

    Return: 

        dict, dict: identities and positives for interactor pairs, respectively

    '''
    seqsim = pd.read_table(InPathSeqSim)
    seqsim['Protein_1'] = seqsim['Protein_1'].map(IDMappingDict)
    seqsim['Protein_2'] = seqsim['Protein_2'].map(IDMappingDict)
    seqsim['Protein_1'], seqsim['Protein_2'] = seqsim[['Protein_1', 'Protein_2']].min(axis=1), seqsim[['Protein_1', 'Protein_2']].max(axis=1)
    seqsim['index'] = seqsim[['Protein_1', 'Protein_2']].agg('_'.join, axis=1)
    identitiesDict = dict([(key, value) for key, value in zip(seqsim['index'], seqsim['Identities'])])
    positivesDict = dict([(key, value) for key, value in zip(seqsim['index'], seqsim['Positives'])])
    return identitiesDict, positivesDict

def import_GI_profile(InPathGI): 
    '''

    This function returns a genetic interaction profile matrix (in DataFrame format). 

    Input: 
        
        InPathGI (Path): Path to genetic interaction profile 

    Return: 

        DataFrame: genetic interaction profile matrix 

    '''
    GI_profile = pd.read_table(InPathGI)
    return GI_profile

def import_coexpression_data(InPathCoexpression): 
    '''

    This function returns gene co-expression profiles (in DataFrame format). 
    If duplicated, the first profile is retained. 

    Input: 
        
        InPathCoexpression (Path): path to co-expression profiles

    Return: 

        DataFrame: gene co-expression profiles

    '''
    cx = pd.read_table(InPath)
    cx = cx.drop(0)
    cx = cx.drop(['NAME', 'GWEIGHT'], axis=1)
    cx = cx.drop_duplicates(subset=['YORF'], keep='first')
    cx = cx.set_index('YORF')
    return cx

def import_semantic_similarity(InPathSemSim, IDMappingFile, index=False): 
    '''

    This function returns semantic similarity (in DataFrame format).

    Input: 

        InPathSemSim (Path): path to semantic similarity

    Return: 

        DataFrame: semantic similarity
    '''
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    semsim = pd.read_table(InPathSemSim)
    semsim.columns = ['Protein_1', 'Protein_2', 'semsim']
    semsim['Protein_1'] = semsim['Protein_1'].map(IDMappingDict)
    semsim['Protein_2'] = semsim['Protein_2'].map(IDMappingDict)
    semsim['Protein_1'], semsim['Protein_2'] = semsim[['Protein_1', 'Protein_2']].min(axis=1), semsim[['Protein_1' ,'Protein_2']].max(axis=1)
    if index: 
        semsim['index'] = semsim[['Protein_1', 'Protein_2']].agg('_'.join, axis=1)
        semsim = semsim.set_index('index')
    return semsim

def change_column_names(InPath, dataset): 
    df = pd.read_table(InPath)
    if dataset == 'HI-Union': 
        df.columns = ['Ensembl_Gene_ID1', 'Ensembl_Gene_ID2', 'Protein_1', 'Protein_2', 'PDB_Structure', 'Protein1_Chain', 'Protein2_Chain', 'Protein1_Interfacial_Residues_Uniprot', 'Protein2_Interfacial_Residues_Uniprot']
    else: 
        df.columns = ['Gene1_Name', 'Gene2_Name', 'Protein1_IntAct_ID', 'Protein2_IntAct_ID', 'Protein_1', 'Protein_2', 'PDB_Structure', 'Protein1_Chain', 'Protein2_Chain', 'Protein1_Interfacial_Residues_Uniprot', 'Protein2_Interfacial_Residues_Uniprot']
    df.to_csv(InPath, sep='\t', index=False)

def get_Ghadie_strInt(InPath, type, OutPath): 
    if type == 'Y2H': 
        df = pd.read_excel(InPath, sheet_name=0)
        df[['Protein_1', 'Protein_2', 'Interfaces']].to_csv(OutPath, sep='\t', index=False)
    elif type == 'Lit': 
        df = pd.read_excel(InPath, sheet_name=2)
        df[['Protein_1', 'Protein_2', 'Interfaces']].to_csv(OutPath, sep='\t', index=False)