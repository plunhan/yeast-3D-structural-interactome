import os
from pathlib import Path
from alphafold_tools import process_alphafold_models

def main(): 

    alphafoldDir = Path('/home/han1204/projects/def-yxia/han1204/alphafold')
    alphafoldModelInfoFile = alphafoldDir / 'alphafold_model_info.xlsx'
    process_alphafold_models(alphafoldDir, alphafoldModelInfoFile)

if __name__ == '__main__':
    main()
