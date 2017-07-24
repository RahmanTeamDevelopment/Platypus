import cython

cimport samtoolsWrapper
cimport fastafile

from samtoolsWrapper cimport Samfile
from samtoolsWrapper cimport IteratorRow
from samtoolsWrapper cimport AlignedRead
from samtoolsWrapper cimport cAlignedRead
from fastafile cimport FastaFile

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

cdef list getHaplotypesInWindow(dict window, int nReads, FastaFile refFile, int ploidy, int maxCoverage, int minMapQual, int minReadQual, int maxHaplotypes, int maxVariants, int maxReadLength)
