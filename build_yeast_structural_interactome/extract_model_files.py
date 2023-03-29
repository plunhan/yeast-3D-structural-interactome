#----------------------------------------------------------------------------------------
# Extract structural model files with extension .ent from a folder.
#----------------------------------------------------------------------------------------

import os
import subprocess
from pathlib import Path

def main():
    
    # input directory containing structural models to extract
    modelDir = Path('../data/processed/HI-II-14/model_based/ppi_models')
    
    # output directory to save extracted models
    outDir = Path('../models_only')
    
    if not modelDir.exists():
        print("Model directory not found")
        return
    
    if not outDir.exists():
        os.makedirs(outDir)
    
    modelFiles = os.listdir(modelDir)
    for name in modelFiles:
        if name.endswith('.ent'):
            cmd = ['cp', str(modelDir / name), str(outDir / name)]
            subprocess.run(cmd)

if __name__ == "__main__":
    main()