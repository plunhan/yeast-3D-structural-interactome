#----------------------------------------------------------------------------------------
# Build an interface-annotated structural interactome from available PPI structural models.
# Interfaces are mapped onto protein sequences.
#----------------------------------------------------------------------------------------

import os
import pandas as pd
import pickle
from pathlib import Path
from modelling_tools import set_model_dir
from pdb_tools import (get_distance, 
                       load_pdbtools_chain_strucRes_labels,
                       load_pdbtools_chain_sequences,
                       structured_chain_residues, 
                       return_chain_res_IDsToPos, 
                       return_structure, 
                       return_chain_res_posToID, 
                       chain_residues)

def structured_residues (pdbid, chainID, pdbDir, resList): 
    '''Return assigned SEQRES residues of a chain. 

    Args: 
        pdbid (str): PDB ID for model parent structure. 
        chainID (str): chain ID. 
        pdbDir (Path): path to file directory containing PDB structures. 
        resList (list): a list of residue IDs on the chain. 

    Returns: 
        dict: list of assigned SEQRES residues

    '''
    struc = return_structure (pdbid, pdbDir)
    if struc: 
        model = struc[0]
        if model.has_id(chainID): 
            residues = chain_residues (model, chainID)
            return [res for res in residues if res.get_id() in resList]
    return []

def calculate_interacting_residues (structuralInteractomeFile, 
                                    pdbDir, 
                                    chainSeqFile, 
                                    chainStrucResFile, 
                                    bindingDist, 
                                    outPath): 
    load_pdbtools_chain_strucRes_labels(chainStrucResFile)
    load_pdbtools_chain_sequences(chainSeqFile)
    print('\t' + 'reading structural interactome')
    interactome = pd.read_table(structuralInteractomeFile)
    interactome['Chain_pairs'] = interactome['Chain_pairs'].apply(lambda x: x.split('|')[0]).apply(lambda x: tuple(map(str.strip, x.split('+'))))
    interactomeWithInterfacialResPairs = []
    for i, ppi in interactome.iterrows(): 
        chain_id1, chain_id2 = ppi.Chain_pairs
        pdbID1, id1 = chain_id1.split('_')
        pdbID2, id2 = chain_id2.split('_')
        pdbID = pdbID1
        interfacialRes1, interfacialRes2 = ppi.Interfaces.split('+')[0].split(','), ppi.Interfaces.split('+')[1].split(',')
        interfacialResPos1 = [return_chain_res_posToID(pdbID, id1, int(resPos), pdbDir) for resPos in interfacialRes1]
        interfacialResPos2 = [return_chain_res_posToID(pdbID, id2, int(resPos), pdbDir) for resPos in interfacialRes2]
        strucInterfacialResPos1 = structured_residues(pdbID1, id1, pdbDir, interfacialResPos1)
        strucInterfacialResPos2 = structured_residues(pdbID2, id2, pdbDir, interfacialResPos2)
        interactingResPairs = [(res1.get_id()[1], res2.get_id()[1]) for res1 in strucInterfacialResPos1 for res2 in strucInterfacialResPos2 if get_distance(res1, res2) < bindingDist]
        interactomeWithInterfacialResPairs.append(interactingResPairs)
        print("\r{:.2f}% has been processed. ".format(len(interactomeWithInterfacialResPairs) * 100 / len(interactome)), end='')
    interactome["Interacting residue pairs"] = interactomeWithInterfacialResPairs
    pickle.dump(interactome, open(outPath, 'wb'))

def main():

    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'Yeast'

    # max binding distance for interface residues in PDB structure
    bindingDist = 5

    # show figures
    showFigs = False

    # parent directory of all data files
    # dataDir = Path('/Volumes/MG_Samsung/edgotype_fitness_effect_full_model/data')
    dataDir = Path('../data')

    # parent directory of all processed data files
    procDir = dataDir / 'processed'

    # directory of aim 2
    aim2Dir = procDir / 'aim_2'

    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name

    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'

    # directory for model structure files
    modelDir = modelBasedDir / 'ppi_models'

    # input data files
    structuralInteractomeFile = modelBasedDir / 'structural_interactome.txt'
    chainSeqFile = modelBasedDir / 'ppi_chain_sequences.pkl'
    modelInterfaceFile = modelBasedDir / 'model_interfaces.txt'
    chainStrucResFile = modelBasedDir / 'ppi_chain_strucRes.pkl'

    # output data files
    interactingResidues = aim2Dir / 'interacting_residues.pkl'
    
    if not modelBasedDir.exists():
        os.makedirs(modelBasedDir)
    
    set_model_dir (modelDir)
    
    if not interactingResidues.is_file(): 
        calculate_interacting_residues (structuralInteractomeFile, 
                                        modelDir, 
                                        chainSeqFile, 
                                        chainStrucResFile, 
                                        bindingDist, 
                                        interactingResidues) 

if __name__ == "__main__":
    main()
