#----------------------------------------------------------------------------------------
# Produce PPI template files required for PPI structural modelling.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from id_mapping import produce_chainSeq_dict
from interactome_tools import write_interactome_sequences
from modelling_tools import (set_pdb_dir,
                             set_template_dir,
                             enable_pdb_downloads,
                             disable_pdb_warnings,
                             write_ppi_template_sequences,
                             produce_fullmodel_chain_strucRes_dict,
                             produce_ppi_template_files)

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct, pombe
    interactome_name = 'pombe'
    
    # allow downloading of PDB structures
    allow_pdb_downloads = False
    
    # suppress PDB warnings
    suppress_pdb_warnings = True
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed template-related data files specific to interactome
    templateBasedDir = interactomeDir / 'template_based'
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # directory for PDB structure files
    pdbDir = Path('../../pdb_files')
    
    # directory for template structure files
    templateDir = modelBasedDir / 'ppi_templates'
    
    # input data files
    # uniqueGeneSequenceFile = procDir / 'yeast_unique_gene_reference_sequences.txt'
    uniqueGeneSequenceFile = procDir / 'pombe_unique_gene_reference_sequences.txt'
    interactomeFile = templateBasedDir / 'structural_interactome.txt'
    chainSeqFile = templateBasedDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = templateBasedDir / 'protein_chain_strucRes.pkl'
    
    # output data files
    interactomeSequenceFile = modelBasedDir / 'protein_sequences.fasta'
    templateSeqFastaFile = modelBasedDir / 'ppi_template_sequences.fasta'
    templateSeqFile = modelBasedDir / 'ppi_template_sequences.pkl'
    templateStrucResFile = modelBasedDir / 'ppi_template_strucRes.pkl'
    
    # create output directories if not existing
    if not modelBasedDir.exists():
        os.makedirs(modelBasedDir)
    if not pdbDir.exists():
        os.makedirs(pdbDir)
    if not templateDir.exists():
        os.makedirs(templateDir)
    
    # set directory of raw PDB coordinate files for modelling tools
    set_pdb_dir (pdbDir)
    
    # set directory of template coordinate files for modelling tools
    set_template_dir (templateDir)
    
    # enable or disable PDB downloads
    enable_pdb_downloads (allow_pdb_downloads)
    
    # suppress or allow PDB warnings
    disable_pdb_warnings (suppress_pdb_warnings)
    
    print('writing PPI protein sequences')
    write_interactome_sequences (interactomeFile,
                                 uniqueGeneSequenceFile,
                                 interactomeSequenceFile)
    
    print('Extracting coordinate files for PPI templates')
    produce_ppi_template_files (interactomeFile, chainSeqFile, chainStrucResFile)
    
    print('writing PPI template sequences to Fasta file')
    write_ppi_template_sequences (interactomeFile,
                                  chainSeqFile,
                                  chainStrucResFile,
                                  templateSeqFastaFile)
    
    print('producing PPI template sequence dictionary')
    produce_chainSeq_dict (templateSeqFastaFile, templateSeqFile)
    
    print('producing PPI template structured residue label file')
    produce_fullmodel_chain_strucRes_dict (templateSeqFile, templateStrucResFile)

if __name__ == "__main__":
    main()
