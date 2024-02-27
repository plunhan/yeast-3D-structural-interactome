import os
import subprocess
import pandas as pd
from pathlib import Path
import pickle
import pymol2
from Bio import PDB
from Bio.PDB import Selection

def process_proteins(entry):
    proteins = entry.split('; ')
    return tuple(proteins)

def map_alphafold_model(alphafoldModelInfoFile): 
    df = pd.read_excel(alphafoldModelInfoFile, usecols=[0,2])
    df['UniProtKB'] = df['UniProtKB'].apply(process_proteins)
    result_dict = df.set_index('ModelArchive ID')['UniProtKB'].to_dict()
    return result_dict

def process_alphafold_models(directory, alphafoldModelInfoFile): 
    alphafoldMappingDict = map_alphafold_model(alphafoldModelInfoFile)
    dir_ls = [d for d in os.listdir(directory) if d.startswith('ma-bak-cepc')]
    gzFile = 'model.cif.gz'
    cifFile = 'model.cif'
    failedItems = set()
    for modelDir in dir_ls: 
        pdbFile = '='.join(alphafoldMappingDict[modelDir])
        decompressCommand = f'gunzip {directory}/{modelDir}/{gzFile}'
        moveCommand = f'mv {directory}/{modelDir}/{cifFile} {directory}/{modelDir}.cif'
        if os.path.exists(f'{directory}/{modelDir}/{gzFile}'): 
            try: 
                subprocess.run(decompressCommand, shell=True, check=True)
                print(f'{gzFile} has been decompressed successfully in {modelDir}. ')
            except subprocess.CalledProcessError as e:
                failedItems.add(modelDir)
                print(f'Error: {e}')
        if os.path.exists(f'{directory}/{modelDir}/{cifFile}'):
            try:
                subprocess.run(moveCommand, shell=True)
                print(f'{gzFile} has been moved to the parent directory. ')
            except subprocess.CalledProcessError as e:
                failedItems.add(modelDir)
                print(f'Error: {e}')
        if os.path.exists(f'{directory}/{modelDir}.cif'): 
            try:
                with pymol2.PyMOL() as p: 
                    p.cmd.load(f'{directory}/{modelDir}.cif', 'myprotein')
                    p.cmd.save(f'{directory}/{pdbFile}.pdb', selection='myprotein')
                print(f'The .cif file has been converted to .pdb file. ')
            except subprocess.CalledProcessError as e:
                failedItems.add(modelDir)
                print(f'Error: {e}')
    print('In total, there are {} models that cannot be converted to .pdb file: '.format(len(failedItems)))
    print('\n'.join(list(failedItems)))

def calculate_interfacial_residues(structure, chain_id1, chain_id2, cutoff_distance=5.0):
    interfacial_residues_1 = set()
    interfacial_residues_2 = set()
    chain1 = structure[0][chain_id1]
    chain2 = structure[0][chain_id2]
    for atom1 in chain1.get_atoms():
        for atom2 in chain2.get_atoms():
            distance = atom1 - atom2
            if distance <= cutoff_distance:
                interfacial_residues_1.add(atom1.parent.id)
                interfacial_residues_2.add(atom2.parent.id)
    if len(interfacial_residues_1) == 0: 
        return [], []
    r1, r2 = sorted([res[1] for res in list(interfacial_residues_1)]), sorted([res[1] for res in list(interfacial_residues_2)])
    r1 = [str(res) for res in r1]
    r2 = [str(res) for res in r2]
    return r1, r2

def get_interfacial_residues(pdbFile, cutoff_distance=5.0):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein_complex', pdbFile)
    chain_id1 = 'A'
    chain_id2 = 'B'
    interfacialRes1, interfacialRes2 = calculate_interfacial_residues(structure, chain_id1, chain_id2, cutoff_distance)
    return interfacialRes1, interfacialRes2

def produce_alphafold_structural_interactome(pdbList, outPath, cutoff_distance=5.0): 
    result_ls = []
    for pdbFile in pdbList: 
        complexName = pdbFile.split('.')[0]
        protein1, protein2 = pdbFile.split('.')[0].split('=')
        res_ls1, res_ls2 = get_interfacial_residues(pdbFile, cutoff_distance)
        if len(res_ls1) == 0: 
            continue
        interfaces = ','.join(res_ls1) + '+' + ','.join(res_ls2)
        result_ls.append( (protein1, protein2, 'Not Applicable', 'AlphaFold2', 'Not Applicable', interfaces, f'{complexName}_A+{complexName}_B') )
        print(f'{pdbFile} has been processed.')
    result_df = pd.DataFrame(result_ls, columns=['Protein_1', 'Protein_2', 'Complex_ID', 'Template_file_ID', 'Alignment_file_ID', 'Interfaces', 'Chain_pairs'])
    result_df.to_csv(outPath, sep='\t')
