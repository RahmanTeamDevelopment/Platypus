#cython: boundscheck=False
#cython: cdivision=True
#cython: nonecheck=False

"""
Fast implementation of the Genotype class.
To profile, put "# cython: profile=True" at the top of this file.
"""

from __future__ import division
import logging
import cython
cimport cython

cimport bamfileutils
cimport chaplotype
cimport samtoolsWrapper

from samtoolsWrapper cimport AlignedRead
from chaplotype cimport Haplotype

###################################################################################################

logger = logging.getLogger("Log")

###################################################################################################

cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten

###################################################################################################

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)

###################################################################################################

cdef class DiploidGenotype(object):
    """
    This class represents a diploid genotype. It stores two
    haplotypes.
    """
    def __init__(self, int ploidy, list haplotypes):
        """
        Constructor. Takes haplotypes.
        """
        if ploidy != 2 or len(haplotypes) != 2:
            raise StandardError, "Diploid genotype must hava ploidy == 2. Supplied ploidy = %s. nHaplotypes = %s" %(ploidy, len(haplotypes))

        self.haplotypes = haplotypes
        self.hap1 = haplotypes[0]
        self.hap2 = haplotypes[1]

    def __str__(self):
        """
        Generate a string representation of the genotype class
        """
        return str(dict(ploidy=2, nHaplotypes=len(self.haplotypes), haplotypes=self.haplotypes))

    def __repr__(self):
        """
        Generate a string representation of the genotype class
        """
        return self.__str__()

    cdef double calculateDataLikelihood(DiploidGenotype self, cAlignedRead** start, cAlignedRead** end, int individualIndex, int nReads, int nIndividuals):
        """
        """
        cdef double likelihood = 0.0
        cdef Py_ssize_t readIndex = 0
        cdef double like1 = 0.0
        cdef double like2 = 0.0
        cdef double* arr1 = self.hap1.alignReads(individualIndex, start, end, nReads, nIndividuals)
        cdef double* arr2 = self.hap2.alignReads(individualIndex, start, end, nReads, nIndividuals)

        for readIndex from 0 <= readIndex < nReads:
            like1 = arr1[readIndex]
            like2 = arr2[readIndex]
            likelihood += log(0.5*(like1 + like2))

        return likelihood

###################################################################################################

cdef list generateAllGenotypesFromHaplotypeList(int ploidy, list haplotypes):
    """
    Return a list of genotypes, corresponding to all allowed combinations of
    haplotypes in the specified region.

    e.g. if nHaplotypes = 2, and n = 2, then the list returned is
    [(0,0), (0,1), (1,1)]

    i.e. we are allowed 2 sets of haplotype 0, 2 sets of haplotype 1, or one of
    each. This would be appropriate for a single diploid individual if we were
    considering 2 haplotypes.
    """
    cdef list genotypes = []
    cdef int nHaplotypes = len(haplotypes)
    cdef int indexOne = 0
    cdef int indexTwo = 0
    cdef Haplotype haplotypeOne
    cdef Haplotype haplotypeTwo

    if ploidy == 1:
        raise StandardError, "Platypus only works for diploid genotypes"

    elif ploidy == 2:

        for indexOne from 0 <= indexOne < nHaplotypes:
            for indexTwo from indexOne <= indexTwo < nHaplotypes:

                haplotypeOne = haplotypes[indexOne]
                haplotypeTwo = haplotypes[indexTwo]

                genotypes.append(DiploidGenotype(ploidy, [haplotypeOne, haplotypeTwo]))

        return genotypes

    else:
        raise StandardError, "Genotype generator can not currently cope with ploidy > 2"

###################################################################################################
