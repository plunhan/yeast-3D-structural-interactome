import os
from pathlib import Path
from alphafold_tools import produce_alphafold_structural_interactome

def main(): 

    alphafoldDir = Path('../alphafold')
    cutoff_distance=5.0

    # output files
    structuralInteractome = alphafoldDir / 'structural_interactome_alphafold.txt'

    if not structuralInteractome.is_file(): 
        pdbList = ['../alphafold/' + file for file in os.listdir(alphafoldDir) if file.endswith('.pdb')]
        produce_alphafold_structural_interactome(pdbList, structuralInteractome, cutoff_distance)

if __name__ == '__main__':
    main()
