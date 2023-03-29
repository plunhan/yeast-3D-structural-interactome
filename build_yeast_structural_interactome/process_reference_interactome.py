#----------------------------------------------------------------------------------------
# Process reference interactome from PPI dataset.
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
from text_tools import (parse_HI_II_14_interactome,
                        parse_HuRI_interactome,
                        parse_IntAct_interactions)
from id_mapping import produce_protein_interaction_dict
from interactome_tools import (write_interactome_sequences,
                               remove_interactions_reported_once,
                               remove_duplicate_PPIs)

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct, pombe
    interactome_name = 'pombe'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    '''
    Yeast_datafile = procDir / 'yeast_reference_interactome.txt'
    uniprotIDmapFile = procDir / 'to_yeast_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_yeast_proteome.list'
    geneMapFile = procDir / 'to_yeast_geneName_map.pkl'
    uniqueGeneSequenceFile = procDir / 'yeast_unique_gene_reference_sequences.txt'
    '''
    Yeast_datafile = procDir / 'pombe_reference_interactome.txt'
    uniprotIDmapFile = procDir / 'to_pombe_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_pombe_proteome.list'
    geneMapFile = procDir / 'to_pombe_geneName_map.pkl'
    uniqueGeneSequenceFile = procDir / 'pombe_unique_gene_reference_sequences.txt'
    
    # output data files
    interactomeFile = interactomeDir / 'reference_interactome.txt'
    interactomeSequenceFile = interactomeDir / 'interactome_sequences.fasta'
    proteinPartnersFile = interactomeDir / 'protein_interaction_partners.pkl'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    interactome = pd.read_table(Yeast_datafile)
    proteins = list(set(interactome[["Protein_1", "Protein_2"]].values.flatten()))
    print('Interactome size without duplicates: %d PPIs, %d proteins' % (len(interactome), len(proteins)))
    
    if not interactomeFile.is_file():
        interactome.to_csv(interactomeFile, index=False, sep='\t')
    
    if not interactomeSequenceFile.is_file():
        print('writing interactome protein sequences')
        write_interactome_sequences (interactomeFile, uniqueGeneSequenceFile, interactomeSequenceFile)
    
    if not proteinPartnersFile.is_file():
        print('producing protein interaction partners dictionary')
        produce_protein_interaction_dict (interactomeFile, proteinPartnersFile)

if __name__ == "__main__":
    main()
