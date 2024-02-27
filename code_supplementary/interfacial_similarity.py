from pathlib import Path
import pandas as pd
from ialign_tools import test_interfacial_similarity

def main(): 

    # cut-off defining hub proteins
    hubCutOff = 2

    # cut-off defining interacting residues (in Angstrom)
    distance = 5

    # reference interactome name
    interactome_name = 'Yeast'

    # parent directory of all data files
    dataDir = Path('../data')

    # parent directory of ialign
    ialignDir = Path('../ialign/bin')

    # parent directory of all external data files
    extDir = dataDir / 'external'

    # parent directory of all processed data files
    procDir = dataDir / 'processed'

    # parent directory of correlation results
    aim1Dir = procDir / 'aim_1'

    # parent directory of round 1 review results
    round1Dir = procDir / 'round_1'

    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'

    # directory of semantic similarity files
    semSimDir = interactomeDir / 'gosim'

    # directory of PPI models
    ppiModelDir = modelBasedDir / 'ppi_models'

    # directory of TM-align files
    TMalignDir = Path('../TMalign/TMalign')

    # input files
    strIntPath = modelBasedDir / 'structural_interactome.txt'
    degreeStrIntPath = aim1Dir / 'degree_count_structural_interactome.txt'
    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'
    GImatrixPath = procDir / 'sga_GI_rel_matrix.txt'
    ialignPath = ialignDir / 'ialign.pl'
    ialignTempFile = ialignDir / 'temp.txt'

    # output files
    outPath = round1Dir / 'confounding_interfacial_similarity_ialign.txt'

    if not outPath.is_file():
        test_interfacial_similarity(strIntPath, 
                                    degreeStrIntPath, 
                                    IDMappingFile, 
                                    GImatrixPath, 
                                    ialignPath, 
                                    ialignTempFile, 
                                    ppiModelDir, 
                                    outPath, 
                                    hub_cutoff=2, 
                                    distance=5, 
                                    method='tm')

if __name__ == '__main__':
    main()