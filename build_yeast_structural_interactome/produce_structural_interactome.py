#----------------------------------------------------------------------------------------
# Build an interface-annotated structural interactome from available PPI structural models.
# Interfaces are mapped onto protein sequences.
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
from modelling_tools import (set_model_dir,
                             produce_model_annotated_interactome,
                             produce_ppi_fullmodel_chainSeq_dict,
                             produce_ppi_fullmodel_pos_mapping,
                             produce_fullmodel_chain_strucRes_dict)
from structural_annotation import (produce_interface_annotated_interactome,
                                   merge_interactome_interface_annotations)

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'pombe'
    
    # max binding distance for interface residues in PDB structure
    bindingDist = 5
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    # dataDir = Path('/Volumes/MG_Samsung/edgotype_fitness_effect_full_model/data')
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # directory for model structure files
    modelDir = modelBasedDir / 'ppi_models'
    
    # input data files
    templateAnnotatedInteractomeFile = modelBasedDir / 'template_annotated_interactome.txt'
    # proteinSeqFile = procDir / 'yeast_reference_sequences.pkl'
    proteinSeqFile = procDir / 'pombe_reference_sequences.pkl'
    
    # output data files
    modelAnnotatedInteractomeFile = modelBasedDir / 'model_annotated_interactome.txt'
    chainSeqFile = modelBasedDir / 'ppi_chain_sequences.pkl'
    chainMapFile = modelBasedDir / 'struc_interactome_chain_map.txt'
    modelInterfaceFile = modelBasedDir / 'model_interfaces.txt'
    chainStrucResFile = modelBasedDir / 'ppi_chain_strucRes.pkl'
    structuralInteractomeFile1 = modelBasedDir / 'structural_interactome_withDuplicates.txt'
    structuralInteractomeFile = modelBasedDir / 'structural_interactome.txt'
    
    if not modelBasedDir.exists():
        os.makedirs(modelBasedDir)
    
    set_model_dir (modelDir)
    
    if not modelAnnotatedInteractomeFile.is_file():
        print('producing model annotated interactome')
        produce_model_annotated_interactome (templateAnnotatedInteractomeFile,
                                             modelAnnotatedInteractomeFile)
    
    if not chainSeqFile.is_file():
        print('producing model chain sequence dictionary')
        produce_ppi_fullmodel_chainSeq_dict (modelAnnotatedInteractomeFile,
                                             proteinSeqFile,
                                             chainSeqFile)
    
    if not chainMapFile.is_file():
        print('producing model chain position mapping file')
        produce_ppi_fullmodel_pos_mapping (modelAnnotatedInteractomeFile,
                                           chainSeqFile,
                                           chainMapFile)
    
    if not chainStrucResFile.is_file():
        print('producing model chain structured residue label file')
        produce_fullmodel_chain_strucRes_dict (chainSeqFile, chainStrucResFile)
    
    if not structuralInteractomeFile1.is_file():
        print('mapping PPI interfaces from structural models')
        produce_interface_annotated_interactome (modelAnnotatedInteractomeFile,
                                                 modelDir,
                                                 chainSeqFile,
                                                 chainMapFile,
                                                 modelInterfaceFile,
                                                 chainStrucResFile,
                                                 1,
                                                 1,
                                                 False,
                                                 0,
                                                 bindingDist,
                                                 structuralInteractomeFile1,
                                                 downloadPDB = False,
                                                 suppressWarnings = False)
        
        print('merging interface annotations for each PPI')
        merge_interactome_interface_annotations (structuralInteractomeFile1,
                                                 structuralInteractomeFile)
    
    structuralInteractome = pd.read_table (structuralInteractomeFile, sep='\t')
    interactomeProteins = list(set(structuralInteractome[["Protein_1", "Protein_2"]].values.flatten()))
    print('\n' + 'Structural interactome:')
    print('%d PPIs' % len(structuralInteractome))
    print('%d proteins' % len(interactomeProteins))
    print()

if __name__ == "__main__":
    main()
