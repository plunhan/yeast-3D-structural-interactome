# Validate pipeline by using PPIs in gold-standard yeast PDB structures. 

import os 
from pathlib import Path
from validation_tools import calculate_positive_rates

def main(): 

    ''' directories '''

    # reference interactome name
    interactome_name = 'Yeast'

    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name

    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'

    rd1Dir = procDir / 'round_1'

    # directory of newly-generated output files in the progress of validation
    validationDir = interactomeDir / 'validation'

    validationModelDir = validationDir / 'model_based'

    # directory for yeast PDB structure files
    pdbYeastDir = procDir / 'yeast_pdb_files'

    ''' input files '''
    GoldStandard = validationModelDir / 'structural_interactome.txt'
    Modelled_5 = modelBasedDir / 'structural_interactome.txt'
    Template = validationDir / 'structural_interactome.txt'
    strucInteractomeChainMapFile = validationDir / 'struc_interactome_chain_map.txt'
    ppiTempAlignmentFile = modelBasedDir / 'ppi_template_extended_alignments.txt'

    ''' output files '''
    disPosRate = rd1Dir / 'distribution_positive_rate_and_seqsim.txt'

    calculate_positive_rates(GoldStandard, 
                             Modelled_5, 
                             Template, 
                             strucInteractomeChainMapFile, 
                             ppiTempAlignmentFile, 
                             disPosRate)

if __name__ == '__main__':
    main()