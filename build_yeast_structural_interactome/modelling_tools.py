import os
import sys
import io
import warnings
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from simple_tools import reverseTuple
from interactome_tools import read_single_interface_annotated_interactome
from pdb_tools import (allow_pdb_downloads,
                       suppress_pdb_warnings,
                       load_pdbtools_chain_sequences,
                       load_pdbtools_chain_strucRes_labels,
                       produce_chain_struc_sequences,
                       pdbfile_id,
                       pdbfile_name,
                       write_partial_structure,
                       get_chain_IDs,
                       structured_residue_IDs,
                       structured_residue_ID)

# directory for PDB structure files
pdbDir = Path('../pdb_files')

# directory for template structure files
templateDir = Path('../templates')

# directory for alignment files
alignmentDir = Path('../alignments')

# directory for model files
modelDir = Path('../models')

def set_pdb_dir (dir):
    
    global pdbDir
    pdbDir = dir
    if not pdbDir.exists():
        os.makedirs(pdbDir)

def set_template_dir (dir):
    
    global templateDir
    templateDir = dir
    if not templateDir.exists():
        os.makedirs(templateDir)

def set_alignment_dir (dir):
    
    global alignmentDir
    alignmentDir = dir
    if not alignmentDir.exists():
        os.makedirs(alignmentDir)

def set_model_dir (dir):
    
    global modelDir
    modelDir = dir
    if not modelDir.exists():
        os.makedirs(modelDir)

def enable_pdb_downloads (download):
    """Set global variable in pdb_tools module to allow PDB downloads.

    Args:
        download (bool): if true, allow PDB downloads.

    """
    allow_pdb_downloads (download)

def disable_pdb_warnings (suppress):
    """Set global variable in pdb_tools module to suppress PDB warnings.

    Args:
        suppress (bool): if true, suppress PDB warnings.

    """
    suppress_pdb_warnings (suppress)

def extend_alignments (inPath, querySeqFile, subjectSeqFile, outPath):
    
    alignments = pd.read_table (inPath, sep="\t")
    with open(querySeqFile, 'rb') as f:
        querySeq = pickle.load(f)
    with open(subjectSeqFile, 'rb') as f:
        subjectSeq = pickle.load(f)
    
    query_alignments, subject_alignments = [], []
    for _, row in alignments.iterrows():
        Qalign, Salign = extend_alignment (querySeq[row.Query],
                                           subjectSeq[row.Subject],
                                           row.Qseq,
                                           row.Sseq,
                                           row.Qstart,
                                           row.Qend,
                                           row.Sstart,
                                           row.Send)
        query_alignments.append(Qalign)
        subject_alignments.append(Salign)
    
    alignments["Qseq"] = query_alignments
    alignments["Sseq"] = subject_alignments
    alignments.to_csv (outPath, index=False, sep='\t')

def extend_alignment (FullQseq,
                      FullSseq,
                      Qseq,
                      Sseq,
                      Qstart,
                      Qend,
                      Sstart,
                      Send):
    
    leftExt, rightExt = Qstart - 1, len(FullQseq) - Qend
    Qalign = FullQseq[:leftExt] + Qseq
    Salign = ('-' * leftExt) + Sseq
    if rightExt > 0:
        Qalign = Qalign + FullQseq[-rightExt:]
        Salign = Salign + ('-' * rightExt)
    
    leftExt, rightExt = Sstart - 1, len(FullSseq) - Send
    Salign = FullSseq[:leftExt] + Salign
    Qalign = ('-' * leftExt) + Qalign
    if rightExt > 0:
        Salign = Salign + FullSseq[-rightExt:]
        Qalign = Qalign + ('-' * rightExt)
    
    return Qalign, Salign

def write_ppi_template_sequences (interactomeFile,
                                  chainSeqFile,
                                  chainStrucResFile,
                                  outPath):
    
    interactome = read_single_interface_annotated_interactome (interactomeFile)
    templateIDs = []
    for templates in interactome["Chain_pairs"].values:
        for chainID1, chainID2 in templates:
            pdbid, c1 = chainID1.split('_')
            _, c2 = chainID2.split('_')
            templateID = '-'.join([pdbid] + sorted([c1, c2]))
            templateIDs.extend([(templateID, c1),(templateID, c2)])
    write_protein_template_sequences (templateIDs, chainSeqFile, chainStrucResFile, outPath)

def write_protein_template_sequences (templateIDs,
                                      chainSeqFile,
                                      chainStrucResFile,
                                      outPath):

    load_pdbtools_chain_sequences (chainSeqFile)
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    produce_chain_struc_sequences (templateIDs, pdbDir, outPath)

def produce_ppi_template_files (inPath, templateSeqFile, templateStrucResFile):
    
    load_pdbtools_chain_sequences (templateSeqFile)
    load_pdbtools_chain_strucRes_labels (templateStrucResFile)
    interactome = read_single_interface_annotated_interactome (inPath)
    n = len(interactome)
    
    templateFiles = os.listdir(templateDir)
    for i, row in interactome.iterrows():
        sys.stdout.write('  PPI %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        for chain1, chain2 in row.Chain_pairs:
            pdbid, chainID1 = chain1.split('_')
            _, chainID2 = chain2.split('_')
            selectChains = sorted([chainID1, chainID2])
            templateID = '-'.join([pdbid] + selectChains)
            filename = pdbfile_name (templateID)
            if filename not in templateFiles:
                resIDs = {c:structured_residue_IDs (pdbid, c, pdbDir) for c in selectChains}
                write_partial_structure (pdbid,
                                         selectChains,
                                         pdbDir,
                                         templateDir / filename,
                                         resIDs = resIDs)
    print()

def produce_protein_template_files (inPath, templateSeqFile, templateStrucResFile):
    
    load_pdbtools_chain_sequences (templateSeqFile)
    load_pdbtools_chain_strucRes_labels (templateStrucResFile)
    templateMap = pd.read_table(inPath, sep='\t')
    templateIDs = templateMap["Subject"].tolist()
    n = len(templateIDs)
    
    templateFiles = os.listdir(templateDir)
    for i, template in enumerate(templateIDs):
        sys.stdout.write('  Template %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        pdbid, chainID = template.split('_')
        templateID = '-'.join([pdbid, chainID])
        filename = pdbfile_name (templateID)
#         if template == '4rer_B':
#             print('1: 4rer_B')
#             print('pdbid = %s' % pdbid)
#             print('chainid = %s' % chainID)
#             print('templateid = %s' % templateID)
#             print('outfile = %s' % str(outFile))
        if filename not in templateFiles:
            resIDs = {chainID : structured_residue_IDs (pdbid, chainID, pdbDir)}
            write_partial_structure (pdbid,
                                     [chainID],
                                     pdbDir,
                                     templateDir / filename,
                                     resIDs = resIDs)
    print()

def produce_ppi_alignment_files (interactomeFile,
                                 alignmentFile,
                                 templateSeqFile,
                                 templateStrucResFile,
                                 outPath):
    
    load_pdbtools_chain_sequences (templateSeqFile)
    load_pdbtools_chain_strucRes_labels (templateStrucResFile)
    interactome = read_single_interface_annotated_interactome (interactomeFile)
    alignments = pd.read_table (alignmentFile, sep="\t")
    n = len(interactome)
    
    allIDs = []
    for i, row in interactome.iterrows():
        sys.stdout.write('  PPI %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        IDs = create_complex_alignment ([row.Protein_1, row.Protein_2],
                                        row.Chain_pairs,
                                        alignments)
        allIDs.append(IDs)
    print()
    interactome["Complex_ID"], interactome["Template_file_ID"], interactome["Alignment_file_ID"] = zip(* allIDs)
    
    missingAlignments = sum(interactome["Alignment_file_ID"] == '-')
    interactome = interactome [interactome["Alignment_file_ID"] != '-']
    print("%d alignment files created" % len(interactome))
    print("%d alignment files not created (template alignment missing)" % missingAlignments)
    interactome.drop(["Interfaces",	"Chain_pairs"], axis=1, inplace=True)
    interactome.to_csv (outPath, index=False, sep='\t')

def produce_protein_alignment_files (templateMapFile,
                                     alignmentFile,
                                     templateSeqFile,
                                     templateStrucResFile,
                                     outPath):
    
    load_pdbtools_chain_sequences (templateSeqFile)
    load_pdbtools_chain_strucRes_labels (templateStrucResFile)
    templateMap = pd.read_table (templateMapFile)
    alignments = pd.read_table (alignmentFile, sep="\t")
    n = len(templateMap)
    
    allIDs = []
    for i, row in templateMap.iterrows():
        sys.stdout.write('  Protein %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        IDs = create_complex_alignment ([row.Query], [(row.Subject,)], alignments)
        allIDs.append(IDs)
    print()
    templateMap["Complex_ID"], templateMap["Template_file_ID"], templateMap["Alignment_file_ID"] = zip(* allIDs)
    missingAlignments = sum(templateMap["Alignment_file_ID"] == '-')
    templateMap = templateMap [templateMap["Alignment_file_ID"] != '-']
    print("%d alignment files created" % len(templateMap))
    print("%d alignment files not created (template alignment missing)" % missingAlignments)
    templateMap.to_csv (outPath, index=False, sep='\t')

def create_complex_alignment (proteins, templates, alignments):
    
    for template in templates:
        pdbid, _ = template[0].split('_')
        chainIDs = [id.split('_')[1] for id in template]
        templateMap = {c:p for c, p in zip(chainIDs, proteins)}
        templateID = '-'.join([pdbid] + sorted(chainIDs))
        
#         chainIDs = get_chain_IDs ('1yxq', Path('/Volumes/MG_Samsung/pdb_files'))
#         print(chainIDs)
#         print('right here')
#         print(get_chain_IDs ('1yxq', Path('/Volumes/MG_Samsung/pdb_files')))
#         print(len(structured_residue_IDs ('1yxq', 'A', Path('/Volumes/MG_Samsung/pdb_files'))))
#         print('right here again')
        chainIDs = get_chain_IDs (templateID, templateDir)
#         print(templateFileID)
#         print(chainIDs)
        orderedProteins = [templateMap[c] for c in chainIDs]
        complexID = '='.join(orderedProteins)
        
        pr_alignments, ch_alignments = [], []
        for protein, chainID in zip(orderedProteins, chainIDs):
            align = alignments [(alignments["Query"] == protein) & 
                                (alignments["Subject"] == '_'.join([templateID, chainID]))]
            if not align.empty:
                prAlign, chAlign = align[["Qseq", "Sseq"]].values[0]
                pr_alignments.append(prAlign)
                ch_alignments.append(chAlign)
#         print(template)
#         print(chainIDs)
        if len(pr_alignments) == len(orderedProteins):
            templateFileID = pdbfile_id (templateID)
            alignmentFileID = '_'.join([complexID, pdbfile_id('-'.join([pdbid] + chainIDs))])
            write_alignment (complexID,
                             templateID,
                             chainIDs,
                             [pr_alignments, ch_alignments],
                             alignmentDir / alignmentFileID)
            return complexID, templateFileID, alignmentFileID
    return '-', '-', '-'

def write_alignment (complexID,
                     templateID,
                     chainIDs,
                     alignments,
                     outPath):

    pr_alignments, ch_alignments = alignments
#     print(templateFileID)
#     print(chainIDs)
#     print(templateDir)
    het, resNum, icode = structured_residue_ID ('first', templateID, chainIDs[0], templateDir)
    firstResID = str(resNum) if icode == ' ' else str(resNum) + icode
    het, resNum, icode = structured_residue_ID ('last', templateID, chainIDs[-1], templateDir)
    lastResID = str(resNum) if icode == ' ' else str(resNum) + icode
    with io.open(outPath, "w") as fout:
        templateFileID = pdbfile_id (templateID)
        fout.write ('>P1;%s\n' % templateFileID)
        fout.write ('structure:%s:%s:%s:%s:%s::::\n' % (templateFileID,
                                                        firstResID,
                                                        chainIDs[0],
                                                        lastResID,
                                                        chainIDs[-1]))
        fout.write ('/'.join(ch_alignments) + '*\n\n')
        fout.write ('>P1;%s\n' % complexID)
        fout.write ('sequence:%s:.:.:.:.::::\n' % complexID)
        fout.write ('/'.join(pr_alignments) + '*\n')

def produce_model_annotated_interactome (inPath, outPath):
    
    allow_pdb_downloads (False)
    suppress_pdb_warnings (False)
    
    interactome = pd.read_table (inPath)
    n, mappingChains = len(interactome), []
    for i, (p1, p2, modelID) in enumerate(interactome[["Protein_1", "Protein_2", "Complex_ID"]].values):
        sys.stdout.write('  PPI %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        chainLetters = get_chain_IDs (modelID, modelDir)
        if chainLetters:
            chainIDs = ['_'.join([modelID, c]) for c in chainLetters]
            complexP1, complexP2 = modelID.split('=')
            if complexP1 == p1:
                mappingChains.append('-'.join(chainIDs))
            else:
                mappingChains.append('-'.join(reverseTuple(chainIDs)))
        else:
            mappingChains.append('-')
    print()
    interactome["Mapping_chains"] = mappingChains
    interactome = interactome [interactome["Mapping_chains"] != '-']
    interactome.to_csv (outPath, index=False, sep='\t')

def produce_ppi_fullmodel_chainSeq_dict (inPath, proteinSeqFile, outPath):
    
    interactome = pd.read_table (inPath)
    with open(proteinSeqFile, 'rb') as f:
        proteinSeq = pickle.load(f)
    
    chainSeq = {}
    for p1, p2, chainMap in interactome[["Protein_1", "Protein_2", "Mapping_chains"]].values:
        c1, c2 = chainMap.split('-')
        for c, p in [(c1, p1), (c2, p2)]:
            if p in proteinSeq:
                chainSeq[c] = proteinSeq[p]
    
    with open(outPath, 'wb') as fOut:
        pickle.dump(chainSeq, fOut)

def produce_protein_fullmodel_chainSeq_dict (inPath, proteinSeqFile, outPath):
    
    templateMap = pd.read_table (inPath)
    with open(proteinSeqFile, 'rb') as f:
        proteinSeq = pickle.load(f)
    
    chainSeq = {}
    for p, m in templateMap[["Query", "Complex_ID"]].values:
        if p in proteinSeq:
            chainSeq[m + '_A'] = proteinSeq[p]
    
    with open(outPath, 'wb') as fOut:
        pickle.dump(chainSeq, fOut)

def produce_ppi_fullmodel_pos_mapping (inPath, chainSeqFile, outPath):
    
    interactome = pd.read_table (inPath)
    with open(chainSeqFile, 'rb') as f:
        chainSeq = pickle.load(f)
    
    mapping = []
    for p1, p2, chainMap in interactome[["Protein_1", "Protein_2", "Mapping_chains"]].values:
        c1, c2 = chainMap.split('-')
        for c, p in [(c1, p1), (c2, p2)]:
            if c in chainSeq:
                pLen = cLen = str(len(chainSeq[c]))
                pPos = cPos = ','.join(map(str, np.arange(1, len(chainSeq[c]) + 1)))
                mapping.append((p, pLen, c, cLen, pPos, cPos))
    
    with io.open(outPath, "w") as fout:
        fout.write ('\t'.join(["Query", "Qlen", "Subject", "Slen", "Qpos", "Spos"]) + '\n')
        for line in mapping:
            fout.write('\t'.join(line) + '\n')

def produce_protein_fullmodel_pos_mapping (inPath, chainSeqFile, outPath):
    
    templateMap = pd.read_table (inPath)
    with open(chainSeqFile, 'rb') as f:
        chainSeq = pickle.load(f)
    
    mapping = []
    for p, m in templateMap[["Query", "Complex_ID"]].values:
        c = m + '_A'
        if c in chainSeq:
            pLen = cLen = str(len(chainSeq[c]))
            pPos = cPos = ','.join(map(str, np.arange(1, len(chainSeq[c]) + 1)))
            mapping.append((p, pLen, c, cLen, '0', '1', '1', pPos, cPos))
    
    with io.open(outPath, "w") as fout:
        fout.write ('\t'.join(["Query", "Qlen", "Subject", "Slen", "Expect", "Qcov", "Scov", "Qpos", "Spos"]) + '\n')
        for line in mapping:
            fout.write('\t'.join(line) + '\n')

def produce_fullmodel_chain_strucRes_dict (chainSeqFile, outPath):
    
    with open(chainSeqFile, 'rb') as f:
        chainSeq = pickle.load(f)
    
    labels = {}
    for id, seq in chainSeq.items():
        labels[id] = '-' * len(seq)
    
    with open(outPath, 'wb') as fOut:
        pickle.dump(labels, fOut)
