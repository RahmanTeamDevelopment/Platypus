import logging

from fastafile cimport FastaFile
from platypus.chaplotype cimport Haplotype
from itertools import combinations


cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)


cdef extern from "math.h":
    double log(double)
    int pow(int, int)


logger = logging.getLogger("Log")


cdef list getFilteredHaplotypes(Haplotype refHaplotype, list variants, int nVars, FastaFile refFile, str windowChr, int windowStart, int windowEnd, int maxHaplotypes, int maxReadLength):
    """
    Now generate the list of all valid haplotypes, given the set of allowed variants,
    and return the top 'maxHaplotypes' valid haplotypes.
    """
    cdef list haplotypeCandidates = []
    cdef tuple varsThisHap
    cdef Haplotype newHaplotype
    cdef double logLikelihoodThisHaplotype
    cdef int nVarsInHap = 0
    cdef int nPossHaplotypes = 0

    # Optimisation: If we only have one variant in the window, the just return the haplotype
    # composed of that variant.
    if nVars == 1:
        return [Haplotype(windowChr, windowStart, windowEnd, (variants[0],), refFile, maxReadLength)]

    # Optimisation: If we only have two variant in the window, the just return all the haplotypes
    # that we can make with these. Check that the combined haplotype is valid.
    elif nVars == 2:
        hap1 = Haplotype(windowChr, windowStart, windowEnd, (variants[0],), refFile, maxReadLength)
        hap2 = Haplotype(windowChr, windowStart, windowEnd, (variants[1],), refFile, maxReadLength)

        if variants[0].refPos != variants[1].refPos:
            hap3 = Haplotype(windowChr, windowStart, windowEnd, (variants[0], variants[1]), refFile, maxReadLength)

            if hap3.valid():
                return [hap1, hap2, hap3]

        return [hap1, hap2]

    # Otherwise we look at all possible combinations, using itertools.combinations.
    else:
        # If there are few enough variants that we can consider all haplotypes, then do
        # that, and return all of them.
        nPossHaplotypes = pow(2, nVars) - 1

        if nPossHaplotypes <= maxHaplotypes:
            for nVarsInHap from 1 <= nVarsInHap <= nVars:

                for varsThisHap in combinations(variants, nVarsInHap):

                    newHaplotype = Haplotype(windowChr, windowStart, windowEnd, varsThisHap, refFile, maxReadLength)

                    if newHaplotype.valid():
                        haplotypeCandidates.append(newHaplotype)
                    else:
                        pass
                        #logger.debug("Created invalid haplotype %s" %(newHaplotype))

            return haplotypeCandidates

        # Otherwise, filter by likelihood.
        else:
            logger.warning("Too many haplotypes (%s when limit is %s). Skipping this window." %(nPossHaplotypes, maxHaplotypes))
            return []
