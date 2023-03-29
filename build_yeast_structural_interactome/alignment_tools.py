#----------------------------------------------------------------------------------------
# Modules for working on sequence alignments.
#----------------------------------------------------------------------------------------

import pickle
import pandas as pd

def extend_alignments (inPath, querySeqFile, subjectSeqFile, outPath):
    """Extend sequence alignments to full length of protein sequences.

    Args:
        inPath (Path): path to tab-delimited file containing sequence alignment table.
        querySeqFile (Path): path to file containing dictionary of query sequences.
        subjectSeqFile (Path): path to file containing dictionary of subject sequences.
        outPath (Path): file path to save extended alignments to.

    """
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
    """Extend a given sequence alignment to full length of protein sequence.

    Args:
        FullQseq (str): full query sequence.
        FullSseq (str): full subject sequence.
        Qseq (str): query side of sequence alignment.
        Sseq (str): subject side of sequence alignment.
        Qstart (int): alignment start position on full query sequence.
        Qend (int): alignment end position on full query sequence.
        Sstart (int): alignment start position on full subject sequence.
        Send (int): alignment end position on full subject sequence.

    """
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
