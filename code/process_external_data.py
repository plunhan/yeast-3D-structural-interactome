# ---------------------------------------------------------------
# Process raw data files from external sources.
# ---------------------------------------------------------------

import os
from pathlib import Path
from data_tools import (process_BIOGRID_protein_interactome, 
                        process_SGA_reliable, 
                        merge_SGA, 
                        generate_GI_matrix, 
                        generate_essential_gene_list, 
                        generate_GI_full_matrix)
from text_tools import parse_fasta_file, reduce_fasta_headers
from id_mapping import (produce_geneName_dict, 
                        produce_uniqueGene_swissProtIDs,
                        produce_uniqueGene_sequences,
                        produce_proteinSeq_dict,
                        produce_uniprotID_dict,
                        produce_rnaToProtein_refseqID_dict)

def main():

    # Parent directory of all data files
    dataDir = Path('../data')
    
    # Parent directory of all external files
    extDir = dataDir / 'external'
    
    # Parent directory of all processed files
    procDir = dataDir / 'processed'

    # Parent directory of SGA pairwise files
    sgaDir = extDir / 'SGA_pairwise'

    # Input files

    # Yeast binary protein-protein interaction datasets
    BIOGRID = extDir / 'BIOGRID-MV-Physical-4.4.200.tab3.txt'

    # Yeast genetic interaction dataset
    sgaEEFile = sgaDir / 'SGA_ExE.txt'
    sgaNNFile = sgaDir / 'SGA_NxN.txt'
    sgaENFile = sgaDir / 'SGA_ExN_NxE.txt'
    sgaDAmPFile = sgaDir / 'SGA_DAmP.txt'

    # ID mapping files
    IDMappingFile = extDir / 'YEAST_559292_idmapping.dat'

    # Uniprot reference sequence file
    uniprotRefSeqFile = extDir / 'UP000002311_559292.fasta'

    # Uniprot reviewed yeast proteome
    proteomeListFile = extDir / 'uniprot_reviewed_yeast_proteome.list'

    # PDB SeqRes file
    pdbSeqresFile = extDir / 'pdb_seqres.txt'

    # Output files

    # Reference binary protein interactome
    proteinInteractome = procDir / 'yeast_reference_interactome.txt'

    # High-quality genetic interactions
    sgaEE_reliable = procDir / 'SGA_ExE_reliable.txt'
    sgaNN_reliable = procDir / 'SGA_NxN_reliable.txt'
    sgaEN_reliable = procDir / 'SGA_ExN_NxE_reliable.txt'
    sgaDAmP_reliable = procDir / 'SGA_DAmP_reliable.txt'
    sgaAll_reliable = procDir / 'SGA_All_reliable.txt'

    # All genetic interactions, regardless of quality
    sgaAll_full = procDir / 'sga_GI_full_matrix.txt'

    # Genetic interaction score matrix, recording genetic interaction scores
    sga_GI_matrix = procDir / 'sga_GI_rel_matrix.txt'

    # Reference sequence fasta file
    refSeqFile = procDir / 'yeast_reference_sequences.fasta'

    # Reference sequence txt file
    sequenceFile = procDir / 'yeast_reference_sequences.txt'

    # ID mapping dictionary (from Uniprot ID to gene name)
    geneMapFile = procDir / 'to_yeast_geneName_map.pkl'

    # Uniprot reviewed yeast proteome
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_yeast_proteome.list'

    # yeast unique gene reference sequences
    uniqueGeneSequenceFile = procDir / 'yeast_unique_gene_reference_sequences.txt'

    # yeast reference protein sequences
    proteinSeqFile = procDir / 'yeast_reference_sequences.pkl'

    # PDB seqres reduced fasta file
    seqresFile = procDir / 'pdb_seqres_reduced.fasta'

    # ID mapping dictionary (from Gene_OrderedLocusName to UniprotAccession)
    uniprotIDmapFile = procDir / 'to_yeast_uniprotID_map.pkl'

    # Essential and non-essential gene lists
    essentialGeneList = procDir / 'essential_genes.pkl'
    nonessentialGeneList = procDir / 'nonessential_genes.pkl'

    if not procDir.exists():
        os.makedirs(procDir)

    if not sgaDir.exists():
        os.makedirs(sgaDir)

    if not proteinInteractome.is_file():
        process_BIOGRID_protein_interactome(BIOGRID, proteinInteractome)

    # Genetic interaction
    if not sgaEE_reliable.is_file():
        print('Processing SGA essential reliable genetic interactions...')
        process_SGA_reliable(sgaEEFile, sgaEE_reliable)
        print('Finished! ')

    if not sgaNN_reliable.is_file():
        print('Processing SGA non-essential reliable genetic interactions...')
        process_SGA_reliable(sgaNNFile, sgaNN_reliable)
        print('Finished! ')

    if not sgaEN_reliable.is_file():
        print('Processing SGA essential vs. non-essential reliable genetic interactions...')
        process_SGA_reliable(sgaENFile, sgaEN_reliable)
        print('Finished! ')

    if not sgaDAmP_reliable.is_file():
        print('Processing SGA DAmP reliable genetic interactions...')
        process_SGA_reliable(sgaDAmPFile, sgaDAmP_reliable)
        print('Finished! ')

    if not sgaAll_reliable.is_file():
        print('Merging SGA reliable genetic interactions into an individual file...')
        merge_SGA([sgaEE_reliable, sgaNN_reliable, sgaEN_reliable], sgaAll_reliable)
        print('Finished! ')

    if not sga_GI_matrix.is_file(): 
        print('Generating SGA reliable genetic interaction matrix...')
        generate_GI_matrix(sgaAll_reliable, sga_GI_matrix)
        print('Finished! ')

    if not sgaAll_full.is_file(): 
        print('Generating SGA full genetic interaction matrix...')
        generate_GI_full_matrix([sgaEEFile, sgaNNFile, sgaENFile], sgaAll_full)
        print('Finished! ')

    # Process files for homologous modelling (adapted from Mohamed Ghadie's Nature Comms 2019)
    # Reduce headers for reference fasta file
    if not refSeqFile.is_file():
        print('Reducing headers in protein sequence fasta file')
        reduce_fasta_headers (uniprotRefSeqFile, '|', 2, 2, refSeqFile)
    # Parse reference fasta file
    if not sequenceFile.is_file():
        print('reading protein sequence fasta file')
        parse_fasta_file (refSeqFile, sequenceFile)
    # Produce Uniprot ID to gene name dictionary
    if not geneMapFile.is_file():
        print('producing UniProtID-to-geneName dictionary')
        produce_geneName_dict (IDMappingFile, proteomeListFile, geneMapFile)
    # Produce unique gene swissProt IDs
    if not uniqueGeneSwissProtIDFile.is_file():
        print('producing list of unique-gene UniProt IDs')
        produce_uniqueGene_swissProtIDs (proteomeListFile, geneMapFile, uniqueGeneSwissProtIDFile)
    # Produce unique gene sequences
    if not uniqueGeneSequenceFile.is_file():
        print('producing sequence file for unique-gene UniProt IDs')
        produce_uniqueGene_sequences (sequenceFile, uniqueGeneSwissProtIDFile, geneMapFile, uniqueGeneSequenceFile)
    
    if not proteinSeqFile.is_file():
        print('producing protein sequence dictionary')
        produce_proteinSeq_dict (refSeqFile, proteinSeqFile)
    
    if not uniprotIDmapFile.is_file():
        print('producing to-UniProt-ID dictionary')
        produce_uniprotID_dict (IDMappingFile, uniqueGeneSwissProtIDFile, uniprotIDmapFile)
    
    if not seqresFile.is_file():
        print('reducing headers in PDB chain sequence file')
        reduce_fasta_headers (pdbSeqresFile, ' ', 1, 1, seqresFile)

    if not essentialGeneList.is_file(): 
        print('generating essential gene list')
        generate_essential_gene_list(essentialGeneList)

if __name__ == '__main__':
    main()