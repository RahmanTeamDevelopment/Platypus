from platypus.samtoolsWrapper cimport cAlignedRead
from platypus.fastafile cimport FastaFile


cdef class bamReadBuffer:
    cdef char* chrom
    cdef int startBase
    cdef int endBase
    cdef int nReads
    cdef int longestRead
    cdef cAlignedRead** cReads
    cdef cAlignedRead** start
    cdef cAlignedRead** end
    cdef cAlignedRead** windowStart
    cdef cAlignedRead** windowEnd
    cdef setWindowPointers(self, int start, int end)


cdef list getHaplotypesInWindow(
        dict window,
        int nReads,
        FastaFile refFile,
        int ploidy,
        int maxCoverage,
        int minMapQual,
        int minReadQual,
        int maxHaplotypes,
        int maxVariants,
        int maxReadLength
)
