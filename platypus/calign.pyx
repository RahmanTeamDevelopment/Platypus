#cython: boundscheck=False
#cython: cdivision=True
#cython: nonecheck=False

"""
Cython module containing glue code for the alignment routines used in
Platypus.
"""

from __future__ import division

cimport cython

import logging
import samtoolsWrapper
cimport cerrormodel
cimport samtoolsWrapper

from samtoolsWrapper cimport cAlignedRead

logger = logging.getLogger("Log")

###################################################################################################

# set the size of the hash.  Ensure that hash_size == math.pow(4,hash_nucs)
cdef int hash_size = 16384
cdef int hash_nucs = 7
cdef int max_sequence_length = hash_size

ctypedef long long size_t
cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten

###################################################################################################

cdef extern from "align.h":
    int fastAlignmentRoutine(char* seq1, char* seq2, char* qual2, int len1, int len2, int gapextend, int nucprior, char* homopolgapq_or_localgapopen, char usehomopolgapq)

###################################################################################################

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *memset(void *buffer, int ch, size_t count )

###################################################################################################

cdef extern from "math.h":
    double exp(double)
    double log(double)

###################################################################################################

#def testAligner():
#    cdef np.ndarray[DTYPE_t] hash
#
#    seq1 = "TATTTGCATGCGCTTTCGAGCTGTTGAAGAGACGTGTATTGGAATAAGTAATCACATAAGTGTTAGTAACTTATTTAAATACGTATAGAGTCGCCTATTTGCCTAGCCTTTTGGTTCTCAGATTTTTTAATTATTACATTGCTATAAGGGTGTAACTGTGTGATAGCCAAAATTTTAAGCTGCAAATGGTTTGTAAATATGATATATTACAAGCTTCATGAAAATCGGTTTATGACTGATCCGCGATTACGTTGAAAGGCGACTGGCAGAGATACTTTTGTTCAGATGTTTTTTCAGGTAGCGATTCCAATGAATAGGTAAAATACCTTGCAAGTTTTGTTGTTGTCGTTGGAGGAAATGTGGATGTGGTTGTTATTGTTGA"
#
#    # set gap open penalties for homopolymers of different lengths; the last penalty is repeated forever
#    homopolgapQ = [50,50,2]
#
#    # set mismatch penalty
#    mismatchQ = 60
#
#    # gap extension penalty
#    gapextend = 5
#
#    # nuc prior
#    nucprior = 5
#
#    # convert penalties to model
#    indel_error_model = {'nonrepetitive': chr(33 + homopolgapQ[0]),
#                         1: ''.join([chr(33+q) for q in homopolgapQ])}
#
#    seq2list = [("AAGTGTTAGTAAC TTATTT AAATAC  GTATAGAGTCGCCTATTTGCC" ,0,0,0),    # 0 mismatches, 0 indels
#                ("AAGTGTTAGTAAC TTATTT AAATAC  GTAcAGAGTCGCCTATTTGCC" ,1,0,0),    # 1 mismatch
#                ("AAGTGTTAGTAAC TTATTT AAATAC TGTATAGAGTCGCCTATTTGCC",0,1,1),     # 1 insertion (homopol 1) of length 1
#                ("AAGTGTTAGTAAC TTATTT AAATAC   TATAGAGTCGCCTATTTGCC"  ,0,1,-1),  # 1 deletion (homopol 1) of length 1
#                ("AAGTGTTAGTAAC TTATTT AAATACTTGTATAGAGTCGCCTATTTGCC",0,1,2),     # 1 insertion (homopol 1) length 2
#                ("AAGTGTTAGTAAC TTATTT AAATAC    ATAGAGTCGCCTATTTGCC"   ,0,1,-2), # 1 deletion (homopol 1) length 2
#                ("AAGTGTTAGTAACTTTATTT AAATAC  GTATAGAGTCGCCTATTTGCC",0,2,1),     # 1 insertion (homopol 2), length 1
#                ("AAGTGTTAGTAAC  TATTT AAATAC  GTATAGAGTCGCCTATTTGCC"  ,0,2,-1),  # 1 deletion (homopol 2), length 1
#                ("AAGTGTTAGTAAC TTATTTAAAATAC  GTATAGAGTCGCCTATTTGCC",0,3,1),     # 1 insertion (homopol 3), length 1
#                ("AAGTGTTAGTAAC TTATTT  AATAC  GTATAGAGTCGCCTATTTGCC"  ,0,3,-1),  # 1 deletion (homopol 3), length 1
#                ]
#
#    s = seq2list[0][0].replace(' ','').upper()
#    #for i in range(1,len(s)):
#    #    seq2list.append( (s[:i]+"AAA"+s[i:], 0,1,-1) )
#
#    for seq2, mismatches, indelcontext, indelsize in seq2list:
#        seq2 = seq2.replace(' ','').upper()
#        len1 = len(seq1)
#        len2 = len(seq2)
#        hash = hash_sequence(seq1)
#        indexOfReadIntoHap = max(0, map_read(seq2, len1, hash) - 8)
#        lenOfHapSeqToTest = len2 + 15
#        qual2 = chr(mismatchQ) * len2
#        homopolgapq = ''.join( [chr(33+q) for q in homopolgapQ] + [chr(0)] )
#        localgapopen = cerrormodel.annotate_sequence( seq1, indel_error_model, 0, 24 )  # output_q_base=0, max_tandem_length=24
#        seq1part = seq1[indexOfReadIntoHap:]
#        localgapopenpart = localgapopen[indexOfReadIntoHap:]
#        score = fastAlignmentRoutine(seq1part, seq2, qual2, lenOfHapSeqToTest, len2, gapextend, nucprior, homopolgapq, True)
#        score2 = fastAlignmentRoutine(seq1part, seq2, qual2, lenOfHapSeqToTest, len2, gapextend, nucprior, localgapopenpart, False)
#        print mismatches, indelcontext, indelsize, score, score2


###################################################################################################

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

###################################################################################################

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

###################################################################################################

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

###################################################################################################

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

###################################################################################################

#cdef void extend_segment( list segment, bytes read, bytes reference, int reference_length ):
#    """
#    Maximally extends segment until first mismatch
#    """
#    cdef int displacement, readstart, refstart, length
#    cdef int readlen = len(read)
#    displacement,readstart,refstart,length = segment
#
#    while (readstart + length < readlen and
#           refstart + length < reference_length and
#           read[readstart+length] == reference[refstart + length]):
#        length += 1
#
#    while (readstart>0 and refstart>0 and
#           read[readstart-1] == reference[refstart-1]):
#        length += 1
#        readstart -= 1
#        refstart -= 1
#        segment[1] -= 1
#        segment[2] -= 1
#    segment[3] = length
#
####################################################################################################
#
#cdef remove_overlap(int a1, int b1, int l1, int a2, int b2, int l2):
#    """
#    Removes any overlap between the alignment blocks, which may have arisen because of their extension
#    """
#    while a1+l1 > a2 or b1+l1 > b2:
#        if l1 > l2:
#            l1 -= 1
#        else:
#            l2 -= 1
#            a2 += 1
#            b2 += 1
#
#    return (a1,b1,l1,a2,b2,l2)
#
####################################################################################################
#
#cpdef map_cortex(bytes read, bytes reference, int reference_length, np.ndarray[DTYPE_t] reference_hash, int min_anchor_length):
#    """
#    Returns an alignment of the (long) cortex read to the reference sequence
#    The alignment is strictly in the forward direction
#    The alignment is encoded as
#       (read_pos_1,ref_pos_1,block_len_1,read_pos_2,ref_pos_2,block_len_2)
#    In case just one aligning segment is identified, the alignment is encoded as
#       (read_pos_1,ref_pos_1,block_len_1)
#    If the read is too short, or no aligning segment is found, None is returned.
#    """
#    assert reference_length < max_sequence_length
#
#    cdef int readlen = len(read)
#    cdef np.ndarray[DTYPE_t] counts = np.zeros([reference_length + readlen], dtype=DTYPE)
#    cdef np.ndarray[DTYPE_t] starts = np.zeros([reference_length + readlen], dtype=DTYPE)
#    cdef char* char_read = read
#    cdef int i, pos, displacement, count, readstart
#    cdef list segments = []
#    cdef int max_gap = 0
#    cdef int best_segment_pair = -1
#    cdef int max_len = -1
#    cdef int best_segment = -1
#
#    if readlen < hash_nucs:
#        return None
#
#    for i in range(readlen - hash_nucs):
#        pos = reference_hash[ hash(char_read + i) ]
#        if pos > 0:
#            displacement = pos - (1+i)
#            assert displacement < reference_length
#            count = counts[<unsigned int>(displacement + readlen)] + 1
#            counts[<unsigned int>(displacement + readlen)] = count
#            if count == 1:
#                # record start of potential new segment
#                starts[<unsigned int>(displacement + readlen)] = i
#            else:
#                if count >= min_anchor_length:
#                    # enter or update segment
#                    readstart = starts[<unsigned int>(displacement+ readlen)]
#                    if len(segments)==0 or segments[-1][0] != displacement:
#                        seqstart = displacement + readstart
#                        segments.append( [displacement,readstart,seqstart,0] )
#                        # if these segments bound the largest gap yet, record
#                        if len(segments)>1 and abs(displacement - segments[-2][0]) > max_gap:
#                            max_gap = abs(displacement - segments[-2][0])
#                            best_segment_pair = len(segments)-2
#                    # update length of segment
#                    segments[-1][3] = i-readstart+1
#    if best_segment_pair == -1:
#        # no large gap found -- return best single segment
#        for i in range(len(segments)):
#            if segments[i][3] > max_len:
#                max_len = segments[i][3]
#                best_segment = i
#        if best_segment == -1:
#            # no segment found
#            return None
#        extend_segment( segments[best_segment], read, reference, reference_length )
#        return ( segments[best_segment][1], segments[best_segment][2], segments[best_segment][3] )
#    # return best pair
#    extend_segment( segments[best_segment_pair], read, reference, reference_length )
#    extend_segment( segments[best_segment_pair+1], read, reference, reference_length )
#    return remove_overlap( segments[best_segment_pair][1], segments[best_segment_pair][2], segments[best_segment_pair][3],
#                           segments[best_segment_pair+1][1], segments[best_segment_pair+1][2], segments[best_segment_pair+1][3] )
#
###################################################################################################
