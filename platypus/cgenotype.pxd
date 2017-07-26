
from chaplotype cimport Haplotype
from platypus.samtoolsWrapper cimport cAlignedRead


cdef class DiploidGenotype:
    cdef public list haplotypes
    cdef Haplotype hap1
    cdef Haplotype hap2
    cdef double calculateDataLikelihood(
            DiploidGenotype self,
            cAlignedRead** start,
            cAlignedRead** end,
            int individualIndex,
            int nReads,
            int nIndividuals
    )


cdef list generateAllGenotypesFromHaplotypeList(int ploidy, list haplotypes)
