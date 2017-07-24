#cython: boundscheck=False
#cython: cdivision=True
#cython: nonecheck=False

"""
Various classes and functions for use in generating and processing
variant candidates.
"""

from __future__ import division

import logging
import fastafile
import cython

cimport chaplotype
cimport variant
cimport fastafile

from variant cimport Variant
from chaplotype cimport convertVariantToStandardFormat
from fastafile cimport FastaFile

logger = logging.getLogger("Log")

###################################################################################################

cdef list filterVariants(list varList, FastaFile refFile, int maxReadLength, int minSupport, int strandFilter, int maxDiff, int verbosity):
    """
    Generator function, this calls the Candidates generator function, and gets a list of
    sorted variant candidates. This list is then merged such that only unique variant candidates
    are returned, and supporting reads are accumulated in the returned candidates.

    No additional filtering is done at this step.
    """
    cdef Variant v
    cdef Variant lastVariant
    cdef Variant newVcfVariant
    cdef list filteredVariants = []
    cdef int support = 0
    cdef int lengthChange = 0

    for v in varList:

        newVcfVariant = convertVariantToStandardFormat(v, refFile, maxReadLength)

        # Start.
        if lastVariant is None:
            lastVariant = newVcfVariant

        # Same as last variant: add supporting reads.
        elif newVcfVariant == lastVariant:
            lastVariant.addVariant(newVcfVariant)

        # Not the same: add to filtered list, if it passes filters.
        else:
            support = lastVariant.nSupportingReads
            lengthChange = abs(len(lastVariant.added) - len(lastVariant.removed))

            if support < minSupport:
                pass
            elif lengthChange > maxDiff:
                pass
            elif strandFilter and (lastVariant.nFwdReads == 0 or lastVariant.nRevReads == 0):
                pass
            else:
                if verbosity >= 3:
                    logger.debug("Adding variant %s to filtered list" %(lastVariant))

                filteredVariants.append(lastVariant)

            lastVariant = newVcfVariant

    if lastVariant is not None:

        support = lastVariant.nSupportingReads
        lengthChange = abs(len(lastVariant.added) - len(lastVariant.removed))

        if support < minSupport:
            pass
        elif lengthChange > maxDiff:
            pass
        elif strandFilter and (lastVariant.nFwdReads == 0 or lastVariant.nRevReads == 0):
            pass
        else:
            if verbosity >= 3:
                logger.debug("Adding variant %s to filtered list" %(lastVariant))

            filteredVariants.append(lastVariant)

    return sorted(filteredVariants)

###################################################################################################
