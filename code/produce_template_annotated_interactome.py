#----------------------------------------------------------------------------------------
# Build a template-annotated interactome from a chain-annotated reference interactome by 
# identifying chain pair annotations that have binding interfaces. Interfaces are also
# mapped onto protein sequences.
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
from pdb_tools import download_structures
from text_tools import parse_blast_file, produce_item_list
from id_mapping import (produce_chainSeq_dict, 
                        produce_chain_strucRes_dict,
                        produce_protein_chain_dict)
from interactome_tools import (read_chain_annotated_interactome,
                               read_single_interface_annotated_interactome)
from structural_annotation import (locate_alignments,
                                   filter_chain_annotations,
                                   produce_alignment_evalue_dict,
                                   produce_chain_annotated_interactome,
                                   produce_interface_annotated_interactome,
                                   merge_interactome_interface_annotations,
                                   remove_duplicate_interface_annotations,
                                   filter_chain_annotations_by_protein)

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct, pombe
    interactome_name = 'Yeast'
    
    # Maximum alignment e-value allowed for protein-chain annotations
    evalue = 1e-5
    
    # Minimum protein coverage fraction required for protein-chain annotation
    proteinCov = 0
    
    # Minimum chain coverage fraction required for protein-chain annotation
    chainCov = 0
    
    # If True, chain residue is required to match aligned protein residue for their 
    # positions to be considered aligned
    resMatch = False
    
    # Minimum fraction of interface residues required to map from chain-pair annotation onto PPI
    mapCutoff = 0.5
    
    # max binding distance for interface residues in template structure
    bindingDist = 5
    
    # If True, merge all mapped interfaces of a PPI into one interface (keep as True)
    mergeInterfaces = True
    
    # max number of templates with mapping interfaces selected for each PPI
    maxTemplates = 2
    
    # max number of chain-pair interface calculations per PPI, including previously known interfaces
    maxAttempts = 100
    
    # If True, randomly sample from a PPI's chain-pair annotations, otherwise start
    # from first annotation (keep as False to select pair with smallest alignment e-value)
    randChainPairs = False
    
    # download missing PDB structures for annotated chain-pairs
    download_missing_structures = True
    
    # allow downloading of PDB structures while constructing the structural interactome
    allow_pdb_downloads = True
    
    # suppress PDB warnings when constructing the structural interactome
    suppress_pdb_warnings = True
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed template-related data files specific to interactome
    templateBasedDir = interactomeDir / 'template_based'
        
    # directory for PDB structure files
    pdbDir = Path('../../pdb_files')
    
    # input data files
    chainResAnnotFile = extDir / 'ss_dis.txt'
    seqresFile = procDir / 'pdb_seqres_reduced.fasta'
    pdbBlastFile = interactomeDir / 'interactome_pdb_e-5'
    interactomeFile = interactomeDir / 'reference_interactome.txt'
    
    # output data files
    chainInterfaceFile = procDir / 'chain_interfaces.txt'
    chainMapFile1 = interactomeDir / 'interactome_chain_alignment.txt'
    chainMapFile2 = interactomeDir / 'interactome_chain_alignment_filtered.txt'
    chainMapFile3 = interactomeDir / 'interactome_chain_map.txt'
    chainSeqFile = templateBasedDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = templateBasedDir / 'protein_chain_strucRes.pkl'
    chainListFile = templateBasedDir / 'protein_model_chains.list'
    modelChainsFile = templateBasedDir / 'protein_model_chains.pkl'
    alignmentEvalueFile = templateBasedDir / 'protein_chain_min_alignment_evalues.pkl'
    chainAnnotatedInteractomeFile = templateBasedDir / 'chain_annotated_interactome.txt'
    chainIDFile = templateBasedDir / 'interactome_chainIDs.txt'
    pdbIDFile = templateBasedDir / 'interactome_pdbIDs.txt'
    templateAnnotatedInteractomeFile1 = templateBasedDir / 'structural_interactome_withDuplicates.txt'
    templateAnnotatedInteractomeFile = templateBasedDir / 'structural_interactome.txt'
    refInteractomeChainMapFile = templateBasedDir / 'ref_interactome_chain_map.txt'
    strucInteractomeChainMapFile = templateBasedDir / 'struc_interactome_chain_map.txt'
    
    if not pdbDir.exists():
        os.makedirs(pdbDir)
    if not templateBasedDir.exists():
        os.makedirs(templateBasedDir)
    
    interactome = pd.read_table(interactomeFile)
    print('\n' + 'Reference interactome:')
    print('%d PPIs' % len(interactome))
    print('%d proteins' % len(set(interactome[["Protein_1", "Protein_2"]].values.flatten())))
    
    if not chainSeqFile.is_file():
        print('producing PDB chain sequence dictionary from fasta records')
        produce_chainSeq_dict (seqresFile, chainSeqFile)
    
    if not chainStrucResFile.is_file():
        print('parsing PDB chain structured-residue order file')
        produce_chain_strucRes_dict (chainResAnnotFile, chainStrucResFile)
    
    if not chainMapFile1.is_file():
        print('parsing BLAST protein-chain alignment file')
        parse_blast_file (pdbBlastFile, chainMapFile1)
    
    if not chainMapFile2.is_file():
        print('filtering chain annotations')
        filter_chain_annotations (chainMapFile1,
                                  chainMapFile2,
                                  evalue = evalue,
                                  prCov = proteinCov,
                                  chCov = chainCov)
    
    if not chainMapFile3.is_file():
        # pausing time in seconds between locating alignments on queries and subjects, occurs only once
        pausetime = 0
        print('locating aligned residues on protein and chain sequences')
        locate_alignments (chainMapFile2,
                           chainMapFile3,
                           resMatch = resMatch,
                           pausetime = pausetime)
    
        print('producing chain ID list')
        produce_item_list (chainMapFile3, "Subject", chainListFile)
    
    if not modelChainsFile.is_file():
        print('producing protein chains dictionary')
        produce_protein_chain_dict (chainMapFile3, modelChainsFile)
    
    if not alignmentEvalueFile.is_file():
        print('producing protein-chain alignment evalue dictionary')
        produce_alignment_evalue_dict (chainMapFile3, alignmentEvalueFile, method = 'min')
    
    if not chainAnnotatedInteractomeFile.is_file():
        print('producing chain-annotated interactome')
        produce_chain_annotated_interactome (interactomeFile,
                                             modelChainsFile,
                                             chainAnnotatedInteractomeFile,
                                             alignmentEvalueFile = alignmentEvalueFile)
        
        chainAnnotatedInteractome = read_chain_annotated_interactome (chainAnnotatedInteractomeFile)
        interactomeProteins = list(set(chainAnnotatedInteractome[["Protein_1", "Protein_2"]].values.flatten()))
        print('Chain-annotated interactome:')
        print('%d PPIs' % len(chainAnnotatedInteractome))
        print('%d proteins' % len(interactomeProteins))
        
        uniqueChains = set()
        for ls in chainAnnotatedInteractome["Mapping_chains"].values:
            for pair in ls:
                uniqueChains.update(pair)
        uniquePDBs = {id.split('_')[0] for id in uniqueChains}
        
        with open(chainIDFile, 'w') as f:
            for i in sorted(uniqueChains):
                f.write("%s\n" % i)
        with open(pdbIDFile, 'w') as f:
            for i in sorted(uniquePDBs):
                f.write("%s\n" % i)
    
    if not refInteractomeChainMapFile.is_file():
        print('filtering chain annotations by reference chain-annotated interactome proteins')
        filter_chain_annotations_by_protein (chainMapFile3,
                                             interactomeProteins,
                                             refInteractomeChainMapFile)
    
    if download_missing_structures:
        print('downloading missing structures for PDB IDs mapping onto interactome')
        download_structures (pdbIDFile, pdbDir)
    
    if not templateAnnotatedInteractomeFile1.is_file():
        print('producing template-annotated interactome')
        produce_interface_annotated_interactome (chainAnnotatedInteractomeFile,
                                                 pdbDir,
                                                 chainSeqFile,
                                                 refInteractomeChainMapFile,
                                                 chainInterfaceFile,
                                                 chainStrucResFile,
                                                 maxTemplates,
                                                 maxAttempts,
                                                 randChainPairs,
                                                 mapCutoff,
                                                 bindingDist,
                                                 templateAnnotatedInteractomeFile1,
                                                 downloadPDB = allow_pdb_downloads,
                                                 suppressWarnings = suppress_pdb_warnings)
        if mergeInterfaces:
            print('merging interface annotations for each PPI')
            merge_interactome_interface_annotations (templateAnnotatedInteractomeFile1,
                                                     templateAnnotatedInteractomeFile)
        else:
            print('removing duplicate interface annotations for each PPI without merging')
            remove_duplicate_interface_annotations (templateAnnotatedInteractomeFile1,
                                                    templateAnnotatedInteractomeFile)
    
    templateInteractome = read_single_interface_annotated_interactome( templateAnnotatedInteractomeFile )
    interactomeProteins = list(set(templateInteractome[["Protein_1", "Protein_2"]].values.flatten()))
    print('\n' + 'Template-annotated interactome:')
    print('%d PPIs' % len(templateInteractome))
    print('%d proteins' % len(interactomeProteins))
    print()
    
    if not strucInteractomeChainMapFile.is_file():
        print('filtering chain annotations by template-annotated interactome proteins')
        filter_chain_annotations_by_protein (chainMapFile3,
                                             interactomeProteins,
                                             strucInteractomeChainMapFile)

if __name__ == "__main__":
    main()
