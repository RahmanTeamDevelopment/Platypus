
from platypus.variant cimport Variant


cdef class Population:
    cdef dict vcfInfo
    cdef dict vcfFilter
    cdef dict vcfVariants
    cdef public list genotypeCalls
    cdef public list genotypeLikelihoods
    cdef public list haplotypes
    cdef double cutoff
    cdef double** genotypePosteriors
    cdef double* frequencies
    cdef int** haplotypeIndexes
    cdef int* nReads
    cdef int nHaplotypes
    cdef int nIndividuals
    cdef int nGenotypes
    cdef int verbosity
    cdef int badReadsWindow
    cdef int badReadsThreshold
    cdef double sbThreshold
    cdef double abThreshold
    cdef list readBuffers

    cdef double calculatePosterior(self, Variant var)
    cdef dict vcfINFO(self)
    cdef dict vcfFILTER(self)
    cdef dict vcfVARIANTS(self)


cdef class Caller:
    cdef int nHaplotypes
    cdef int nIndividuals
    cdef int nGenotypes
    cdef int ploidy
    cdef int verbosity
    cdef list haplotypes
    cdef list genotypes
    cdef list readBuffers
    cdef double** genotypePosteriors
    cdef double** cisr
    cdef int** haplotypeIndexes
    cdef int* nReads
    cdef object options
    cdef Population call(self, int maxIters=*)
    cdef double EMiteration(self, double* freq)
    cdef setup(self, list haplotypes, list genotypes, int nInd, int ploidy, int verbosity, list readBuffers)
