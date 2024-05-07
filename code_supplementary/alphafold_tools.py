import os
import subprocess
import pandas as pd
from pathlib import Path
import pickle
import pymol2
from Bio import PDB
from Bio.PDB import Selection
import multiprocessing as mp

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
    residues1 = chain1.get_residues()
    residues2 = chain2.get_residues()
    chSeq1 = ''
    for residue in residues1:
        if PDB.is_aa(residue):
            chSeq1 += PDB.Polypeptide.three_to_one(residue.get_resname())
    chSeq2 = ''
    for residue in residues2:
        if PDB.is_aa(residue):
            chSeq2 += PDB.Polypeptide.three_to_one(residue.get_resname())
    for atom1 in chain1.get_atoms():
        for atom2 in chain2.get_atoms():
            distance = atom1 - atom2
            if distance <= cutoff_distance:
                interfacial_residues_1.add(atom1.parent.id)
                interfacial_residues_2.add(atom2.parent.id)
    if len(interfacial_residues_1) == 0: 
        return [], [], chSeq1, chSeq2
    r1, r2 = sorted([res[1] for res in list(interfacial_residues_1)]), sorted([res[1] for res in list(interfacial_residues_2)])
    r1 = [str(res) for res in r1]
    r2 = [str(res) for res in r2]
    return r1, r2, chSeq1, chSeq2

def get_interfacial_residues(pdbFile, cutoff_distance=5.0):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein_complex', pdbFile)
    chain_id1 = 'A'
    chain_id2 = 'B'
    intRes1, intRes2, chSeq1, chSeq2 = calculate_interfacial_residues(structure, chain_id1, chain_id2, cutoff_distance)
    return intRes1, intRes2, chSeq1, chSeq2

def worker(args):
    pdbFile, protSeq_dict, cutoff_distance = args
    pdbPath = pdbFile
    pdbFile = pdbFile.split('/')[-1]
    complexName = pdbFile.split('.')[0]
    protein1, protein2 = pdbFile.split('.')[0].split('=')
    if (protein1 not in protSeq_dict) or (protein2 not in protSeq_dict):
        return ('0', '0', '0', '0', '0', '0', '0')
    protSeq1, protSeq2 = protSeq_dict[protein1], protSeq_dict[protein2]
    res_ls1, res_ls2, chSeq1, chSeq2 = get_interfacial_residues(pdbPath, cutoff_distance)
    if len(res_ls1) == 0: 
        return ('0', '0', '0', '0', '0', '0', '0')
    if (chSeq1 == protSeq2) and (chSeq2 == protSeq1): 
        res_ls1, res_ls2 = res_ls2, res_ls1
    elif (chSeq1 == protSeq1) and (chSeq2 == protSeq2):
        res_ls1, res_ls2 = res_ls1, res_ls2
    else:
        return ('0', '0', '0', '0', '0', '0', '0')
    interfaces = ','.join(res_ls1) + '+' + ','.join(res_ls2)
    return (protein1, protein2, 'Not Applicable', 'AlphaFold2', 'Not Applicable', interfaces, f'{complexName}_A+{complexName}_B')

def produce_alphafold_structural_interactome(pdbList, protSeq_dict, outPath, cutoff_distance=5.0): 
    ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
    pool = mp.Pool(processes=ncpus)
    worker_args = [(pdbFile, protSeq_dict, cutoff_distance) for pdbFile in pdbList]
    result_ls = pool.map(worker, worker_args)
    result_df = pd.DataFrame(result_ls, columns=['Protein_1', 'Protein_2', 'Complex_ID', 'Template_file_ID', 'Alignment_file_ID', 'Interfaces', 'Chain_pairs'])
    result_df.to_csv(outPath, sep='\t')
