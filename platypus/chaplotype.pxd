import cython

cimport fastafile
cimport bamfileutils
cimport samtoolsWrapper
cimport variant

from samtoolsWrapper cimport cAlignedRead
from fastafile cimport FastaFile
from variant cimport Variant

cdef class Haplotype:
    cdef bytes refName
    cdef int startPos
    cdef int endPos
    cdef tuple variants
    cdef bytes haplotypeSequence
    cdef char* cHaplotypeSequence
    cdef FastaFile refFile
    cdef double** likesByIndCache
    cdef int lenCache
    cdef bytes referenceSequence
    cdef int hash
    cdef int hapLen
    cdef int varSeqLen
    cdef char* cHomopolQ
    cdef dict indelErrorModel
    cdef int useIndelErrorModel
    cdef bytes localGapOpenQ
    cdef char* cLocalGapOpenQ
    cdef double* alignReads(self, int individualIndex, cAlignedRead** start, cAlignedRead** end, int nReads, int nIndividuals)
    cdef int valid(self)
    cdef char* getReferenceSequence(self, prefix=*)
    cdef char* getMutatedSequence(self)
    cdef double strandBiasPValue(self, Variant variant)
    cdef dict vcfFILTER(self)
    cdef dict vcfINFO(self)
    cdef bytes getSequenceContext(self, Variant variant)
    cdef int homopolymerLengthForOneVariant(self, Variant variant)
    cdef list homopolymerLengths(self)

cdef Variant convertVariantToStandardFormat(Variant variant, FastaFile refFile, int maxReadLength)
cdef double binomial(int x, int size, double prob)
cdef double logFactorial(int x)
