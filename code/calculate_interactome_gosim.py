#----------------------------------------------------------------------------------------
# Compare functional similarity of interaction partners between the structural 
# interactome, the reference interactome and random interactions. Similarity is 
# calculated by the fastsemsim library using a selected semantic similarity measure
# applied to gene ontology (GO) terms of interaction partners.
#
# Run the following scripts before running this script:
# - process_interactome.py
# - produce_structural_interactome.py
#
# Requirements:
# - install fastsemsim using the command line 'pip install fastsemsim'
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from protein_function import produce_fastsemsim_protein_gosim_dict
import itertools


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
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'Yeast'
    
    # similarity measure to calculate GO similarity
    # options: Resnik, Lin, Jiang-Conrath, SimGIC, SimUI, SimIC, SimRel, Dice, SimTO
    #           SimNTO, Jaccard, Czekanowski-Dice, Cosine, GSESAME, SimICND, SimICNP
    sim_measure = 'SimGIC'
    
    # mixing strategy for merging GO term semantic similarities
    # options: max, avg, BMA (best match average)
    mix_method = 'BMA'
    
    # root ontology on which GO similarity is calculated
    # options: biological_process, molecular_function, cellular_component
    ont_root = 'cellular_component'
    
    # root ontology labels used to label output files and figures
    ont_abv = {'biological_process':'bp', 'molecular_function':'mf', 'cellular_component':'cc'}
    
    # list of ontological relationships to ignore
    ont_ignore = None
    
    # list of evidence codes to ignore
    ec_ignore = None
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of GO similarity output data files
    gosimDir = interactomeDir / 'gosim'
    
    # get root ontology label to label output files and figures
    ont_label = ont_abv [ont_root]
    
    # input data files
    ontologyFile = extDir / 'go-basic.obo'
    annotationFile = extDir / 'goa_yeast.gaf'
    interactomeFile = interactomeDir / 'reference_interactome.txt'
    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'
    
    # output data files
    refPPIListFile = gosimDir / 'refPPIs.txt'
    refPPIgosimParamFile = gosimDir / ('mf_fastsemsim_parameters_refPPIs_%s_%s' % (ont_label, sim_measure))
    refPPIfastsemsimOutFile = gosimDir / ('mf_fastsemsim_output_refPPIs_%s_%s' % (ont_label, sim_measure))
    refPPIgosimFile = gosimDir / ('mf_gosim_refPPIs_%s_%s.pkl' % (ont_label, sim_measure))
    
    # create directories if not existing
    if not gosimDir.exists():
        os.makedirs(gosimDir)

    # create ID mapping dictionary to convert the names of reference PPIs
    IDMappingDict = map_protein_id_to_locus_id(IDMappingFile)
    
    #------------------------------------------------------------------------------------
    # load reference and structural interactomes
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (interactomeFile)
    # interactome = interactome[interactome['Protein_1'].isin(IDMappingDict.keys()) & interactome['Protein_2'].isin(IDMappingDict.keys())]
    # interactome['Protein_1'] = interactome['Protein_1'].map(IDMappingDict)
    # interactome['Protein_2'] = interactome['Protein_2'].map(IDMappingDict)
    interactomeProteins = list(set(interactome[["Protein_1", "Protein_2"]].values.flatten()))
    print( '\n' + 'reference interactome:' )
    print( '%d PPIs' % len(interactome) )
    print( '%d proteins' % len(interactomeProteins) )
        
    #------------------------------------------------------------------------------------
    # Calculate GO similarity for all interaction partners in the reference interactome
    #------------------------------------------------------------------------------------
    
    refPPIs = pd.DataFrame()
    refPPIs["Protein_1"], refPPIs["Protein_2"] = zip(* interactome[["Protein_1","Protein_2"]].values)
    refPPIList = itertools.combinations(list(set(refPPIs['Protein_1']).union(set(refPPIs['Protein_2']))), 2)
    refPPIs = pd.DataFrame(refPPIList, columns=['Protein_1', 'Protein_2'])
    
    # produce protein GO similarity dictionary
    if not refPPIgosimFile.is_file():
        print('\n' + 'producing GO similarity dictionary for every possible pairs between proteins in reference PPIs')
        refPPIs[["Protein_1","Protein_2"]].to_csv(refPPIListFile, index=False, header=False, sep='\t')
        produce_fastsemsim_protein_gosim_dict (refPPIListFile,
                                               refPPIgosimFile,
                                               sim_measure = sim_measure,
                                               mix_method = mix_method,
                                               ont_root = ont_root,
                                               ont_ignore = ont_ignore,
                                               ec_ignore = ec_ignore,
                                               ontologyFile = ontologyFile,
                                               annotationFile = annotationFile,
                                               paramOutFile = refPPIgosimParamFile,
                                               fastsemsimOutFile = refPPIfastsemsimOutFile)
    with open(refPPIgosimFile, 'rb') as f:
        gosim = pickle.load(f)
        
    sim = []
    for p in map(tuple, map(sorted, refPPIs[["Protein_1","Protein_2"]].values)):
        sim.append(gosim[p] if p in gosim else np.nan)
    refPPIs["gosim"] = sim
    refPPIs = refPPIs [np.isnan(refPPIs["gosim"]) == False].reset_index(drop=True)
    
if __name__ == '__main__':
    main()
