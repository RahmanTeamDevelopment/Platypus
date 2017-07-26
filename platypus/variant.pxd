
from platypus.fastafile cimport FastaFile
from platypus.samtoolsWrapper cimport AlignedRead
from platypus.samtoolsWrapper cimport cAlignedRead


cdef class Variant:
    cdef public bytes refName
    cdef public bytes added
    cdef public bytes removed
    cdef public int refPos
    cdef public int minPos
    cdef public int maxPos
    cdef public int nSupportingReads
    cdef public int nFwdReads
    cdef public int nRevReads
    cdef public int hashValue
    cdef public int nAdded
    cdef public int nRemoved
    cdef public int sumSqReadPos
    cdef public double rmsReadPos
    cdef public double nFwdReadsWeighted
    cdef public double nRevReadsWeighted
    cdef public double freqPrior
    cdef double calculatePrior(self)
    cdef void addVariant(self, Variant other)


cdef class VariantCandidateGenerator:
    cdef int CIGAR_M
    cdef int CIGAR_I
    cdef int CIGAR_D
    cdef int CIGAR_N
    cdef int CIGAR_S
    cdef int CIGAR_H
    cdef int CIGAR_P
    cdef int minMapQual
    cdef int minBaseQual
    cdef int minFlank
    cdef int refId
    cdef int maxCoverage
    cdef int verbosity
    cdef int genSNPs
    cdef int genIndels
    cdef int maxReadLength
    cdef long int rStart
    cdef long int rEnd
    cdef long int refSeqStart
    cdef long int refSeqEnd
    cdef char* refSeq
    cdef dict variantHeap
    cdef list bamFiles
    cdef bytes pyRefSeq
    cdef bytes rname
    cdef FastaFile refFile

    cdef addVariantToList(self, Variant var)
    cdef getSnpCandidatesFromReadSegment(self, cAlignedRead* read, char* readSeq, char* readQual, int readStart, int readOffset, int refOffset, int lenSeqToCheck)
    cdef getVariantCandidatesFromSingleRead(self, cAlignedRead* read)
    cdef addCandidatesFromReads(self, cAlignedRead** readStart, cAlignedRead** readEnd)
    cdef list getCandidates(self)

    cpdef getListOfSnpPositions(self, AlignedRead read)
    cpdef getListOfIndelPositions(self, AlignedRead read)
