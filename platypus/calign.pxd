cimport platypus.samtoolsWrapper
from platypus.samtoolsWrapper cimport cAlignedRead


cdef int* hash_sequence(char* sequence)

cdef double alignNoTraceback(
        cAlignedRead* read,
        char* hapSeq,
        int hapStart,
        int varSeqLen,
        char* homopolgapq_or_localgapopen,
        int gapextend,
        int nucprior,
        int hapLen,
        int useAsHomopolGapQ
)
