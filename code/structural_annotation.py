#----------------------------------------------------------------------------------------
# Modules for protein structural annotation.
#----------------------------------------------------------------------------------------

import os
import io
import sys
import time
import warnings
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from random import sample
from simple_tools import create_dir, isolate_pairs, merge_list_pairs
from text_tools import read_list_table
from interactome_tools import (write_chain_annotated_interactome,
                               read_chain_annotated_interactome,
                               write_unmerged_interface_annotated_interactome,
                               read_unmerged_interface_annotated_interactome,
                               write_single_interface_annotated_interactome,
                               write_interface_annotated_interactome)
from pdb_tools import (allow_pdb_downloads,
                       suppress_pdb_warnings,
                       clear_structures,
                       load_pdbtools_chain_strucRes_labels,
                       load_pdbtools_chain_sequences,
                       valid_strucRes,
                       get_interface_by_chainIDs)

#-----------------------------------------
# Global variables modified by modules
#-----------------------------------------

# dictionary of chain interfaces already calculated
known_interfaces = {}

#-----------------------------------------
# Modules
#-----------------------------------------

def load_dictionaries (chainSequenceFile = None, chainStrucResLabelFile = None):
    """Load dictionaries containing PDB chain sequence data.

    Args:
        chainSequenceFile (Path): path to file containing PDB chain sequence dictionary.
        chainStrucResLabelFile (Path): path to file containing PDB chain labels for residues 
                                        with available 3D coordinates.

    """
    if chainSequenceFile:
        load_pdbtools_chain_sequences (chainSequenceFile)
    if chainStrucResLabelFile:
        load_pdbtools_chain_strucRes_labels (chainStrucResLabelFile)

def clear_interfaces ():
    """Clear dict of known chain interfaces
    
    """
    global known_interfaces
    print('\t' + 'clearing known interfaces dictionary')
    known_interfaces.clear()

def produce_alignment_evalue_dict (inPath, outPath, method = 'min'):
    """Produce dictionary of protein-chain sequence alignment e-values.

    Args:
        inPath (Path): path to tab-delimited file containing sequence alignments.
        outPath (Path): file path to save alignment e-values to.
        method (str): method for combining e-values from multiple alignments for the same 
                        protein-chain pair, 'min' or 'mean'.

    """
    chainMap = pd.read_table(inPath, sep='\t')
    evalueDict = {}
    if method is 'min':
        for i, row in chainMap.iterrows():
            k = row.Query + '-' + row.Subject
            if k in evalueDict:
                if (row.Expect < evalueDict[k]):
                    evalueDict[k] = row.Expect
            else:
                evalueDict[k] = row.Expect
    elif method is 'mean':
        for i, row in chainMap.iterrows():
            k = row.Query + '-' + row.Subject
            if k in evalueDict:
                evalueDict[k].append( row.Expect )
            else:
                evalueDict[k] = [ row.Expect ]
        for k in evalueDict:
            evalueDict[k] = np.mean( evalueDict[k] )
    with open(outPath, 'wb') as fOut:
        pickle.dump(evalueDict, fOut)
    
def locate_alignments (inPath,
                       outPath,
                       resMatch = False,
                       pausetime = 0):
    """Transform sequence alignments to sequence residue positions.

    Args:
        inPath (Path): path to tab-delimited file containing sequence alignment table.
        outPath (Path): file path to save processed alignment table to.
        resMatch (bool): if True, residues must match for positions to be considered aligned.
        pausetime (numeric): pausing time in seconds between locating alignments on 
                                protein sequences and chain sequences.

    """
    alignments = pd.read_table(inPath, sep="\t")
    print('\t' + '%d alignments to locate' % len(alignments))
    
    print('\t' + 'locating alignments on query sequences')
    alignments["Qpos"] = alignments.apply(lambda x: locate_alignment(x["Qseq"],
                                                                     x["Sseq"],
                                                                     x["Qstart"],
                                                                     resMatch = resMatch), axis=1)
    
    print('\t' + 'pausing for %d seconds before locating alignments on chain sequences' % pausetime)
    time.sleep(pausetime)
    
    print('\t' + 'locating alignments on PDB sequences')
    alignments["Spos"] = alignments.apply(lambda x: locate_alignment(x["Sseq"],
                                                                     x["Qseq"],
                                                                     x["Sstart"],
                                                                     resMatch = resMatch), axis=1)
    
    alignments.drop(["Qseq","Sseq","Match"], axis=1, inplace=True)
    alignments["Qpos"] = alignments["Qpos"].apply(lambda x: ','.join(map(str, x)))
    alignments["Spos"] = alignments["Spos"].apply(lambda x: ','.join(map(str, x)))
    alignments.to_csv(outPath, index=False, sep='\t')

def locate_alignment (Qseq, Sseq, Qstart, resMatch = False):
    """Transform single sequence alignment to positions on query sequence.

    Args:
        Qseq (str): query sequence.
        Sseq (str): subject sequence.
        Qstart (int): alignment start position on query sequence.
        resMatch (bool): if True, residues must match for positions to be considered aligned.
    
    Returns:
        list

    """
    if resMatch:
        matchPos = [i for i, ch in enumerate(Qseq) if (ch != '-') and (Sseq[i] == ch)]
    else:
        matchPos = [i for i, ch in enumerate(Qseq) if (ch != '-') and (Sseq[i] != '-')]
    gapPos = [i for i, ch in enumerate(Qseq) if ch == '-']
    if len(gapPos) == 0:
        return [pos + Qstart for pos in matchPos]
    else:
        numGaps = [len([g for g in gapPos if g < pos]) for pos in matchPos]
        return [pos + Qstart - gaps for pos, gaps in zip(matchPos, numGaps)]
                                 
def map_positions (refPos, alignPos, pos):
    """Map a list of reference positions to their aligned positions.

    Args:
        refPos (list): reference positions.
        alignPos (list): positions aligning with reference positions.
        pos (list): reference positions to map.
    
    Returns:
        list, numeric: mapped positions, fraction of positions mapped.

    """
    mappedPos = [alignPos[refPos.index(i)] for i in pos if i in refPos]
    frac = round(len(mappedPos) / len(pos), 2)
    return mappedPos, frac

def produce_chain_annotated_interactome (inPath,
                                         proteinChainsFile,
                                         outPath,
                                         alignmentEvalueFile = None):
    """Annotate protein-protein interactions (PPIs) with chain pairs aligned to protein sequences.

    Args:
        inPath (Path): path to tab-deleimited file containing PPIs.
        proteinChainsFile (Path): path to file containing dict of chains aligned to each protein.
        outPath (Path): file path to save chain-annotated interactome.
        alignmentEvalueFile (Path): dictionary of protein-chain sequence alignment e-values.

    """
    annotatedInteractome = pd.read_table(inPath, sep="\t")
    with open(proteinChainsFile, 'rb') as f:
        proteinChains = pickle.load(f)
    if alignmentEvalueFile:
        with open(alignmentEvalueFile, 'rb') as f:
            alignmentEvalues = pickle.load(f)
    else:
        alignmentEvalues = None
    
    annotatedInteractome["Mapping_chains"] = annotatedInteractome.apply(lambda x:
                                                    sameStruct_mapping_chains(x["Protein_1"],
                                                                              x["Protein_2"],
                                                                              proteinChains,
                                                                              alignmentEvalues = alignmentEvalues),
                                                    axis=1)
    
    annotatedInteractome = annotatedInteractome[annotatedInteractome["Mapping_chains"].apply(len)
                                                > 0].reset_index(drop=True)
    print('\t%d interactions annotated with same-structure chains' % len(annotatedInteractome))
    write_chain_annotated_interactome (annotatedInteractome, outPath)
    
def sameStruct_mapping_chains (protein1,
                               protein2,
                               proteinChains,
                               alignmentEvalues = None):
    """Identify same-structure chain pairs mapping onto two interacting proteins.

    Args:
        protein1 (str): ID of first interacting protein.
        protein2 (str): ID of second interacting protein.
        proteinChains (dict): dict of chains mapping to each protein.
        alignmentEvalues (dict): e-values for protein-chain sequence alignment.

    Returns:
        list

    """
    chainPairs = []
    if (protein1 in proteinChains) and (protein2 in proteinChains):
        p1chains = sorted(set(proteinChains[protein1]))
        p2chains = sorted(set(proteinChains[protein2]))
    
        pdbID1 = [c.split('_')[0] for c in p1chains]
        pdbID2 = [c.split('_')[0] for c in p2chains]
        
        for i, c1 in enumerate(p1chains):
            for j, c2 in enumerate(p2chains):
                if (c1 != c2) and (pdbID1[i] == pdbID2[j]):
                    chainPairs.append( (c1, c2) )
        
        if alignmentEvalues:
            avgEvalues = []
            for c1, c2 in chainPairs:
                evalue1 = alignmentEvalues[protein1 + '-' + c1]
                evalue2 = alignmentEvalues[protein2 + '-' + c2]
                avgEvalues.append( np.mean( [evalue1, evalue2] ) )
            sortedPairs = [ chainPairs[i] for i, _ in sorted(enumerate(avgEvalues), key=lambda x: x[1]) ]
            return sortedPairs   
    return chainPairs

def filter_chain_annotations (inPath,
                              outPath,
                              evalue = 1e-5,
                              prCov = 0,
                              chCov = 0):
    """Filter protein chain annotations by coverage and e-value.

    Args:
        inPath (Path): path to tab-deleimited file containing protein-chain alignments.
        outPath (Path): file path to save filtered alignments to.
        evalue (numeric): maximum alignment e-value cutoff.
        prCov (numeric): minimum alignment coverage fraction required on protein sequence.
        chCov (numeric): minimum alignment coverage fraction required on chain sequence.

    """
    with io.open(inPath, "r", encoding="utf8", errors='ignore') as f, io.open(outPath, "w") as fout:
        headers = f.readline().strip()
        fout.write(headers + '\t' + 'Qcov' + '\t' + 'Scov' + '\n')
        headerSplit = headers.split('\t')
        numCol = len(headerSplit)
        evalueCol = headerSplit.index('Expect')
        QseqCol = headerSplit.index('Qseq')
        SseqCol = headerSplit.index('Sseq')
        QlenCol = headerSplit.index('Qlen')
        SlenCol = headerSplit.index('Slen')
    with io.open(inPath, "r", encoding="utf8", errors='ignore') as f, io.open(outPath, "a") as fout:
        next(f)
        for line in f:
            line = line.strip()
            linesplit = line.split('\t')
            if len(linesplit) == numCol:
                if float(linesplit[evalueCol]) < evalue:
                    Qseq = linesplit[QseqCol]
                    Sseq = linesplit[SseqCol]
                    Qcov = (len(Qseq) - Qseq.count('-')) / int(linesplit[QlenCol])
                    Scov = (len(Sseq) - Sseq.count('-')) / int(linesplit[SlenCol])
                    if (Qcov >= prCov) and (Scov >= chCov):
                        fout.write('\t'.join([line, str(Qcov), str(Scov)]) + '\n')
    remove_duplicate_chain_annotations (outPath, outPath)

def remove_duplicate_chain_annotations (inPath, outPath):
    """Keep only one alignment for each protein-chain pair, the one with the smallest 
        e-value.

    Args:
        inPath (Path): path to tab-deleimited file containing protein-chain alignments.
        outPath (Path): file path to save filtered alignments to.

    """
    chainMap = pd.read_table(inPath, sep='\t')
    chainMap = chainMap.sort_values("Expect", axis=0, ascending=True)
    chainMap = chainMap.drop_duplicates(subset = ["Query", "Subject"], keep='first')
    chainMap = chainMap.sort_values(['Query','Subject'], axis=0, ascending=True)
    chainMap.to_csv(outPath, index=False, sep='\t')

def filter_chain_annotations_by_protein (inPath, proteins, outPath):
    """Filter protein chain annotations by proteins.

    Args:
        inPath (Path): path to tab-deleimited file containing protein-chain alignments.
        proteins (list): proteins to be selected.
        outPath (Path): file path to save filtered alignments to.

    """
    with io.open(inPath, "r", errors='ignore') as f, io.open(outPath, "w") as fout:
        headers = f.readline()
        fout.write(headers)
        headerSplit = headers.strip().split('\t')
        numCol = len(headerSplit)
        queryPos = headerSplit.index('Query')
    with io.open(inPath, "r", errors='ignore') as f, io.open(outPath, "a") as fout:
        next(f)
        for line in f:
            linesplit = line.strip().split('\t')
            if len(linesplit) == numCol:
                if linesplit[queryPos] in proteins:
                    fout.write(line)

def single_chain_per_protein (inPath,
                              outPath,
                              chainSeqFile = None,
                              chainStrucResFile = None,
                              pdbDir = None):
    """Keep only one chain alignment for each protein, the one with the smallest e-value,
        and valid structured residue sequence if chain structured resdidue label file is
        provided.

    Args:
        inPath (Path): path to tab-deleimited file containing protein-chain alignments.
        outPath (Path): file path to save filtered alignments to.
        chainSeqFile (Path): path to file containing dictionary of model chain sequences.
        chainStrucResFile (Path): path to file containing dict of labels for chain sequence 
                                    positions associated with 3D coordinates.
        pdbDir (Path): file directory containing PDB structures.

    """
    chainMap = pd.read_table (inPath, sep='\t')
    chainMap = chainMap.sort_values ("Expect", axis=0, ascending=True)
    if chainSeqFile and chainStrucResFile and pdbDir:
        load_pdbtools_chain_sequences (chainSeqFile)
        load_pdbtools_chain_strucRes_labels (chainStrucResFile)
        keep = pd.Series(data = False, index = chainMap.index)
        proteins = list(set(chainMap["Query"].values))
        n = len(proteins)
        for k, p in enumerate(proteins):
            sys.stdout.write('  Protein %d out of %d (%.2f%%) \r' % (k+1, n, 100*(k+1)/n))
            sys.stdout.flush()
            for i, c in chainMap.loc[chainMap["Query"]==p, "Subject"].iteritems():
                pdbid, chainID = c.split('_')
                if valid_strucRes (pdbid, chainID, pdbDir):
                    keep[i] = True
                    break
        chainMap = chainMap[keep].reset_index(drop=True)
    chainMap = chainMap.drop_duplicates (subset = "Query", keep='first')
    chainMap = chainMap.sort_values ("Query", axis=0, ascending=True)
    chainMap.to_csv (outPath, index=False, sep='\t')

def produce_interface_annotated_interactome (inPath,
                                             pdbDir,
                                             chainSeqFile,
                                             chainMapFile,
                                             interfaceFile,
                                             chainStrucResFile,
                                             maxInterfaces,
                                             maxAttempts,
                                             rnd,
                                             mapCutoff,
                                             bindingDist,
                                             outPath,
                                             downloadPDB = True,
                                             suppressWarnings = False):
    """Map interaction interfaces onto chain-annotated interactome.

    Args:
        inPath (Path): path to file containing chain-annotated interactome.
        pdbDir (Path): file directory containing pdb structures.
        chainSeqFile (Path): path to file containing dictionary of model chain sequences.
        chainMapFile (Path): path to tab-delimited file containing protein-chain alignments.
        interfaceFile (Path): path to file containing already calculated chain-pair interfaces.
                                Newly calculated interfaces are added to this file.
        chainStrucResFile (Path): path to file containing dict of labels for chain sequence 
                                    positions associated with 3D coordinates.
        maxInterfaces (int): maximum number of interfaces mapped onto single PPI from structural models.
        maxAttempts (int): maximum number of chain-pair annotations per PPI scanned for interface.
        rnd (bool): if True, randomly select chain pairs from PPI chain-pair annotation list,
                    otherwise done by order.
        mapCutoff (numeric): minimum fraction of interface residues mapping onto protein 
                            required for interface to be selected.
        bindingDist (numeric): maximum distance in Angstroms allowed between interface residues
                                of two binding chains.
        outPath (Path): file path to save interface-annotated interactome to.
        downloadPDB (bool): if True, PDB structure downloads are allowed.
        suppressWarnings (bool): if True, PDB warnings are suppressed.

    """
    global known_interfaces
    clear_structures()
    clear_interfaces()
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    load_dictionaries (chainSequenceFile=chainSeqFile, chainStrucResLabelFile=chainStrucResFile)
    print('\t' + 'reading chain-annotated interactome')
    interactome = read_chain_annotated_interactome(inPath)
    print('\t' + 'reading protein-chain mapping table')
    chainMap = read_list_table(chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    annotatedPPIs, numPPIs = 0, len(interactome)
    
    if interfaceFile.is_file():
        print('\t' + 'loading known chain interfaces')
        read_chain_interfaces(interfaceFile)
    else:
        print('\t' + 'creating interface file directory')
        filedir, _ = os.path.split(interfaceFile)
        create_dir(filedir)
    
    print('\t' + 'resolving interfaces for all %d PPIs' % numPPIs)
    interfaces = []
    for i, ppi in interactome.iterrows():
        print('\t' + 'PPI: ' + ppi.Protein_1 + ' - ' + ppi.Protein_2)
        sel_interfaces, sel_frac, sel_chainpairs = map_multi_interfaces (pdbDir,
                                                                         ppi.Protein_1,
                                                                         ppi.Protein_2,
                                                                         ppi.Mapping_chains,
                                                                         chainMap,
                                                                         maxInterfaces,
                                                                         maxAttempts,
                                                                         rnd,
                                                                         mapCutoff,
                                                                         bindingDist,
                                                                         interfaceFile,
                                                                         PPIcount = i + 1,
                                                                         numPPIs = numPPIs)
        interfaces.append( (sel_interfaces, sel_frac, sel_chainpairs) )
        if sel_chainpairs:
            annotatedPPIs += 1
        print('\t' + '%d PPIs resolved' % (i + 1))
        print('\t' + '%d PPIs annotated' % annotatedPPIs)
    interactome["Interfaces"] = interfaces
    
    print('\t' + 'finalizing interface-annotated interactome')
    interactome = interactome[interactome["Interfaces"].apply(lambda x: len(x[0]) > 0)].reset_index(drop=True)
    interactome.drop("Mapping_chains", axis=1, inplace=True)
    interactome["Chain_pairs"] = interactome["Interfaces"].apply(lambda x: x[2])
    interactome["Mapping_frac"] = interactome["Interfaces"].apply(lambda x: x[1])
    interactome["Interfaces"] = interactome["Interfaces"].apply(lambda x: x[0])
    write_unmerged_interface_annotated_interactome (interactome, outPath)

def map_multi_interfaces (pdbDir,
                          protein1,
                          protein2,
                          chainPairs,
                          chainMap,
                          maxInterfaces,
                          maxAttempts,
                          rnd,
                          mapCutoff,
                          bindingDist,
                          interfaceFile,
                          PPIcount = None,
                          numPPIs = None):
    """Map interaction interfaces onto a protein-protein interaction.

    Args:
        pdbDir (Path): file directory containing pdb structures.
        protein1 (str): UniProt ID of first interacting protein.
        protein2 (str): UniProt ID of second interacting protein.
        chainPairs (list): chain pairs mapping onto the two interacting proteins.
        chainMap (DataFrame): protein-chain alignment table.
        maxInterfaces (int): maximum number of interfaces mapped onto single PPI from structural models.
        maxAttempts (int): maximum number of chain-pair annotations per PPI scanned for interface.
        rnd (bool): if True, randomly select chain pairs from PPI chain-pair annotation list,
                    otherwise done by order.
        mapCutoff (numeric): minimum fraction of interface residues mapping onto protein 
                            required for interface to be selected.
        bindingDist (numeric): maximum distance in Angstroms allowed between interface residues
                                of two binding chains.
        interfaceFile (Path): path to file containing already calculated chain-pair interfaces.
                                Newly calculated interfaces are added to this file.
        PPIcount (int): number of PPI currently being annotated.
        numPPIs (int): total number of PPIs annotated so far.

    Returns:
        list, list, list: mapped interfaces, interface mapping coverage, structural models used.

    """
    numAttempts, interfacesFound, numChainPairs = 0, 0, len(chainPairs)
    sel_interfaces, sel_frac, sel_chainpairs, sel_pdbs = [], [], [], []
    if rnd:
        chPairs = sample(chainPairs, numChainPairs)
    else:
        chPairs = chainPairs
    for i, chainPair in enumerate(chPairs):
        if (numAttempts == maxAttempts) or (interfacesFound == maxInterfaces):
            break
        chainID1, chainID2 = chainPair
        pdbID, _ = chainID1.split('_')
        if PPIcount and numPPIs:
            print('\t\t' + 'PPI: %s - %s (# %d out of %d) (%.1f %%)' 
                    % (protein1, protein2, PPIcount, numPPIs, 100 * PPIcount / numPPIs))
        else:
            print('\t\t' + 'PPI: %s - %s' % (protein1, protein2))
        print('\t\t' + '%d attempts (max %d)' % (numAttempts, maxAttempts))
        print('\t\t' + '%d interfaces mapped (seeking %d)' % (interfacesFound, maxInterfaces))
        print('\t\t' + 'Chain pair: %s - %s (# %d out of %d)' % (chainID1, chainID2, i + 1, numChainPairs))
        if pdbID not in sel_pdbs:
            interfaces1, interfaces2 = map_twoside_interfaces (pdbDir,
                                                               protein1,
                                                               protein2,
                                                               chainPair,
                                                               bindingDist,
                                                               chainMap,
                                                               interfaceFile)
            passFrac1 = [f for _, f in interfaces1 if f >= mapCutoff]
            passFrac2 = [f for _, f in interfaces2 if f >= mapCutoff]
            if (len(passFrac1) > 0) and (len(passFrac2) > 0):
                inf1 = [i for i, f in interfaces1 if f >= mapCutoff]
                inf2 = [i for i, f in interfaces2 if f >= mapCutoff]
                sel_interfaces.append((inf1, inf2))
                sel_frac.append((passFrac1, passFrac2))
                sel_chainpairs.append(([chainID1]*len(inf1), [chainID2]*len(inf2)))
                sel_pdbs.append(pdbID)
                interfacesFound += 1
                print('\t\t\t' + 'interface %d found' % interfacesFound)
            numAttempts += 1
        else:
            print('\t\t\t' + 'interface already mapped from PDB structure')
    return sel_interfaces, sel_frac, sel_chainpairs

def map_twoside_interfaces (pdbDir,
                            protein1,
                            protein2,
                            chainPair,
                            bindingDist,
                            chainMap,
                            interfaceFile):
    """Map interaction interfaces from one chain-pair onto a protein-protein interaction.

    Args:
        pdbDir (Path): file directory containing pdb structures.
        protein1 (str): UniProt ID of first interacting protein.
        protein2 (str): UniProt ID of second interacting protein.
        chainPair (tuple): chain pair mapping onto the two interacting proteins.
        bindingDist (numeric): maximum distance in Angstroms allowed between interface residues
                                of two binding chains.
        chainMap (DataFrame): protein-chain alignment table.
        interfaceFile (Path): path to file containing already calculated chain-pair interfaces.
                                Newly calculated interfaces are added to this file.

    Returns:
        list, list: interfaces on first protein, interfaces on second protein.

    """
    global known_interfaces
    chain1, chain2 = chainPair
    chainKey1 = chain1 + '-' + chain2
    chainKey2 = chain2 + '-' + chain1
    if (chainKey1 not in known_interfaces) or (chainKey2 not in known_interfaces):
        print('\t\t\tcomputing interface for chain pair ' + chainKey1)
        known_interfaces[chainKey1], known_interfaces[chainKey2] = get_interface_by_chainIDs (pdbDir,
                                                                                              chain1,
                                                                                              chain2,
                                                                                              maxDist = bindingDist)
        print('\t\t\t' + 'writing chain interfaces to file')
        write_chain_interfaces(interfaceFile)
    else:
        print('\t\t\tinterface known for chain pair ' + chainKey1)
    if (len(known_interfaces[chainKey1]) > 0) and (len(known_interfaces[chainKey2]) > 0):
        alignment1 = chainMap[(chainMap["Query"]==protein1) &
                              (chainMap["Subject"]==chain1)]
        alignment2 = chainMap[(chainMap["Query"]==protein2) &
                              (chainMap["Subject"]==chain2)]
        interfaces1 = alignment1.apply(lambda x:
                                             map_positions(x["Spos"],
                                                           x["Qpos"],
                                                           known_interfaces[chainKey1]),
                                              axis=1)
        interfaces2 = alignment2.apply(lambda x:
                                             map_positions(x["Spos"],
                                                           x["Qpos"],
                                                           known_interfaces[chainKey2]),
                                              axis=1)
        interfaces1 = [(i, f) for i, f in interfaces1.values if f > 0]
        interfaces2 = [(i, f) for i, f in interfaces2.values if f > 0]
        return interfaces1, interfaces2
    else:
        return [], []

def merge_interactome_interface_annotations (inPath, outPath):
    """Merge interface annotations for each PPI in an interactome into one interface.

    Args:
        inPath (Path): path to file containing interface-annotated interactome.
        outPath (Path): file path to save interactome with merged interfaces to.
    
    """
    interactome = read_unmerged_interface_annotated_interactome ( inPath )
    interactome["Interfaces"] = interactome["Interfaces"].apply( isolate_pairs )
    interactome["Interfaces"] = interactome["Interfaces"].apply( merge_list_pairs )
    interactome["Chain_pairs"] = interactome["Chain_pairs"].apply( isolate_pairs )
    interactome = interactome.drop("Mapping_frac", axis=1)
    write_single_interface_annotated_interactome (interactome, outPath)

def remove_duplicate_interface_annotations (inPath, outPath):
    """Remove duplicate interface annotations per PPI from interactome.

    Args:
        inPath (Path): path to file containing interface-annotated interactome.
        outPath (Path): file path to save processed interface-annotated interactome to.
    
    """
    interactome = read_unmerged_interface_annotated_interactome ( inPath )
    interactome["Interfaces"] = interactome["Interfaces"].apply( isolate_pairs )
    interactome["Mapping_frac"] = interactome["Mapping_frac"].apply( isolate_pairs )
    
    uniqueInterfaces = interactome.apply(lambda x:
                                           drop_duplicate_interfaces(x["Interfaces"],
                                                                     x["Mapping_frac"]),
                                           axis=1)
    interactome["Interfaces"] = uniqueInterfaces.apply(lambda x: x[0])
    interactome = interactome.drop("Mapping_frac", axis=1)
    write_interface_annotated_interactome( outPath )

def drop_duplicate_interfaces (interfaces, fractions):
    """Remove duplicate interfaces from a list of interfaces.

    Args:
        interfaces (list): list of interfaces.
        fractions (list): mapping coverage associated with interfaces.
    
    Returns:
        list, list: unique interfaces, mapping coverages.

    """
    uniqueInt = []
    uniqueFrac = []
    for i, intfc in enumerate(interfaces[:-1]):
        if all([intfc != intfc2 for intfc2 in interfaces[i+1:]]):
            uniqueInt.append(intfc)
            uniqueFrac.append(fractions[i])
    uniqueInt.append(interfaces[-1])
    uniqueFrac.append(fractions[-1])
    return uniqueInt, uniqueFrac

def filter_interfaces_by_frac (interfaces, fractions, cutoff):
    """Filter interfaces from a list of interfaces by mapping coverage.

    Args:
        interfaces (list): list of interfaces.
        fractions (list): mapping coverage associated with interfaces.
        cutoff (numeric): cutoff used for both sides of interface.
    
    Returns:
        list, list: interfaces, mapping coverage.
    
    """
    selInt = [interfaces[i] for i, (a, b) in enumerate(fractions) if (a >= cutoff) and (b >= cutoff)]
    selFrac = [(a, b) for (a, b) in fractions if (a >= cutoff) and (b >= cutoff)]
    return selInt, selFrac

def write_chain_interfaces (outPath):
    """Write chain interfaces to file.

    Args:
        outPath (Path): file path to save chain interfaces to.
    
    """
    with io.open(outPath, "w") as fout:
        fout.write('\t'.join(["Chain_1", "Chain_2", "Chain1_interface"]) + '\n')
        for chainKey in sorted(known_interfaces.keys()):
            chain1, chain2 = chainKey.split('-')
            if len(known_interfaces[chainKey]) > 0:
                interface = ','.join([str(elm) for elm in known_interfaces[chainKey]])
            else:
                interface = '-1'
            fout.write('\t'.join([chain1, chain2, interface]) + '\n')

def read_chain_interfaces (inPath):
    """Read chain interfaces from file.

    Args:
        inPath (Path): path to file containing chain interfaces.
    
    """
    if inPath.is_file():
        known_interfaces_df = read_list_table(inPath, "Chain1_interface", int, '\t')
        for _, row in known_interfaces_df.iterrows():
            chainKey = row.Chain_1 + '-' + row.Chain_2
            if -1 in row.Chain1_interface:
                known_interfaces[chainKey] = []
            else:
                known_interfaces[chainKey] = row.Chain1_interface
