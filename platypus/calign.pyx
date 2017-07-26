from __future__ import division

import logging

cimport platypus.samtoolsWrapper
from platypus.samtoolsWrapper cimport cAlignedRead

logger = logging.getLogger("Log")


# set the size of the hash.  Ensure that hash_size == math.pow(4,hash_nucs)
cdef int hash_size = 16384
cdef int hash_nucs = 7
cdef int max_sequence_length = hash_size

ctypedef long long size_t
cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten


cdef extern from "align.h":
    int fastAlignmentRoutine(char* seq1, char* seq2, char* qual2, int len1, int len2, int gapextend, int nucprior, char* homopolgapq_or_localgapopen, char usehomopolgapq)


cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *memset(void *buffer, int ch, size_t count )


cdef extern from "math.h":
    double exp(double)
    double log(double)


cdef double alignNoTraceback(cAlignedRead* read, char* hapSeq, int hapStart, int varSeqLen, char* homopolgapq_or_localgapopen, int gapextend, int nucprior, int hapLen, int useAsHomopolGapQ):
    """
    This is the basic, banded-alignment routine that forms the heart of Platypus. This function decides where to anchor the
    read sequence to the specified haplotype, and calls the fastAlignmentRoutine function, which performs a banded alignment.

    If we don't anchor the read sequence to the correct part of the haplotype, then all the results, particularly for indels,
    will be rubbish.
    """
    cdef char* readSeq = read.seq
    cdef char* readQuals = read.qual

    cdef int readStart = read.pos
    cdef int readLen = read.rlen
    cdef int mapQual = read.mapq

    cdef int lenOfHapSeqToTest = readLen + 15
    cdef int hapPos = readStart - hapStart
    cdef int alignScore = 0

    cdef double probMapWrong = exp(mLTOT*mapQual)
    cdef double probMapRight = 1.0 - probMapWrong

    # If there are no indels in the haplotype (or total insertion + deletion length == 0), just use the left alignment
    # point as an anchor.
    if varSeqLen == 0:
        assert lenOfHapSeqToTest + max(0, hapPos-8) <= hapLen
        alignScore = fastAlignmentRoutine(hapSeq+max(0, hapPos-8), readSeq, readQuals, lenOfHapSeqToTest, readLen, gapextend, nucprior, homopolgapq_or_localgapopen, useAsHomopolGapQ)
        return min(1.0, ((exp(mLTOT*alignScore) * probMapRight) + probMapWrong) )

    # Otherwise, we'll need to re-map the read to the haplotype
    cdef int* hapHash = hash_sequence(hapSeq)
    cdef int indexOfReadIntoHap = max(0, map_read(readSeq, hapLen, hapHash) - 8)
    cdef bytes temp

    free(hapHash)

    # TODO: Put this somewhere else.
    if not (lenOfHapSeqToTest + indexOfReadIntoHap <= hapLen):
        temp = bytes(hapSeq) + "N"*readLen
        hapSeq = temp

    alignScore = fastAlignmentRoutine(hapSeq+indexOfReadIntoHap, readSeq, readQuals, lenOfHapSeqToTest, readLen, gapextend, nucprior, homopolgapq_or_localgapopen, useAsHomopolGapQ)
    return min(1.0, ((exp(mLTOT*alignScore) * probMapRight) + probMapWrong) )


cdef unsigned int my_hash(char* seq):
    """
    Encodes first hash_nucs nucleotides (A,C,G,T or a,c,g,t) into an integer
    """
    cdef int i, c
    cdef unsigned int h = 0

    for i in range(hash_nucs):
        c = seq[i] & 7   # a,A->1  c,C->3  g,G->7  t,T->4

        if c == 7:
            c = 2

        h = (h << 2) + <unsigned int>(c & 3)

    return h


cdef int* hash_sequence(char* sequence):
    """
    Returns a hash of the byte sequence
    """
    cdef int seqlen = len(sequence)
    cdef char* seq = sequence
    cdef int* h = <int*>(calloc(hash_size, sizeof(int)))
    cdef int i
    cdef unsigned int hidx

    if seqlen < hash_nucs:
        return h

    assert seqlen < max_sequence_length

    for i in range(seqlen - hash_nucs):

        hidx = my_hash(seq+i)

        if h[hidx] != 0:
            h[hidx] = -1
        else:
            h[hidx] = i+1

    return h


cdef map_read(char* read, int sequence_length, int* sequence_hash):
    """
    Map a single read to a small-ish (~1kb) sequence. This is used to find the best anchor point
    in a haplotype sequence for the specified read.
    Returns the read's index into the sequence.  May be negative.  Assumes a forward read direction
    """
    assert sequence_length < max_sequence_length

    cdef int readlen = len(read)
    cdef int* counts = <int*>(calloc(sequence_length + readlen, sizeof(int)))
    cdef int maxcount = 0
    cdef int maxpos = 0
    cdef int i, pos, count

    if readlen < hash_nucs:
        return 0

    for i in range(readlen - hash_nucs):
        pos = sequence_hash[ my_hash(read + i) ]

        if pos > 0:
            pos -= 1+i   # account for +1 at indexing; account for index into read
            assert pos < sequence_length

            count = counts[<unsigned int>(pos + readlen)] + 1
            counts[<unsigned int>(pos + readlen)] = count

            if count > maxcount:
                maxcount = count
                maxpos = pos

    free(counts)
    return maxpos
