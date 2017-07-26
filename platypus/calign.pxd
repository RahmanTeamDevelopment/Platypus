cimport cython
cimport bamfileutils
cimport samtoolsWrapper
from samtoolsWrapper cimport cAlignedRead

#cpdef map_cortex(bytes read, bytes reference, int reference_length, int* reference_hash, int min_anchor_length)

cdef int* hash_sequence(char* sequence)
cdef double alignNoTraceback(cAlignedRead* read, char* hapSeq, int hapStart, int varSeqLen, char* homopolgapq_or_localgapopen, int gapextend, int nucprior, int hapLen, int useAsHomopolGapQ)
