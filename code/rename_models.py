#----------------------------------------------------------------------------------------
# Rename structural models created by modeller as "pdb<modelID>.ent".
#----------------------------------------------------------------------------------------

import os
from pathlib import Path

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct, Yeast
    interactome_name = 'Yeast'
    
    # type of models: ppi_models or protein_models
    modelType = 'ppi_models'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # directory for output models
    modelDir = modelBasedDir / modelType
    
    modelFiles = os.listdir(modelDir)
    for name in modelFiles:
        if name.endswith('.B99990001.pdb'):
            modelID, _, _ = name.split('.')
            os.rename(modelDir / name, modelDir / ('pdb' + modelID + '.ent'))

if __name__ == "__main__":
    main()
