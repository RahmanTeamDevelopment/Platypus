from __future__ import division
import logging


logger = logging.getLogger("Log")


class WindowGenerator(object):
    """
    This class is used to generate a list of all possible variants for a given region.
    """
    def __init__(self, refFile, sourceFile=None):
        """
        Constructor. Takes a bam-file reader, a fasta file of the reference,
        and a window size, which tells it how large a region to consider, for
        each set of haplotypes.
        """
        self.refFile = refFile
        self.sourceFile = sourceFile

    def getVariantsByPos(self, chromosome, start, end, sortedVariants):
        """
        Returns a list of lists of variant candidates. Each list contains all candidates
        at that position.
        """
        varsByPos = {}

        for v in sortedVariants:
            if v.refName == chromosome and v.refPos >= start and v.refPos <= end:
                try:
                    varsByPos[v.refPos].append(v)
                except KeyError:
                    varsByPos[v.refPos] = [v]

        listOfVarsByPos = []

        for pos in sorted(varsByPos.keys()):
            listOfVarsByPos.append(varsByPos[pos])

        return listOfVarsByPos

    def getBunchesOfInteractingVariants(self, varsByPos, readLength):
        """
        Go through the list of lists or variants sorted by position, and
        concatenate neighbouring lists, if any of the variants in those
        lists interact.
        """
        bunchedVars = []

        for varList in varsByPos:

            if len(bunchedVars) == 0:
                bunchedVars.append(varList)
            else:
                maxLastPos = max([var.maxPos for var in bunchedVars[-1]])
                minThisPos = min([var.minPos for var in varList])

                # Always put interacting variants in the same window
                if maxLastPos >= minThisPos:
                    bunchedVars[-1].extend(varList)
                # If there is no interaction, start a new window
                else:
                    bunchedVars.append(varList)

        return bunchedVars

    def getWindowVariants(self, chromosome, start, end, sortedVariants, readLength, maxVar):
        """
        Iterate over the list of bunches of interacting vars, and concatenate
        neighbouring bunches, if the total number of variants in the two neighbouring
        bunches is <= maxVarsInWindow. Return a list of lists, where each sub-list
        is a set of variants to be considered in the same window.
        """
        varsByPos = self.getVariantsByPos(chromosome, start, end, sortedVariants)
        bunchedVars = self.getBunchesOfInteractingVariants(varsByPos, readLength)
        return bunchedVars

    def WindowsAndVariants(self, chromosome, start, end, sortedVariants, options):
        """
        Generator: return a series of dictionaries, each containing the start and end
        positions of a haplotype window, plus a list of all the variants in that window.
        Basically takes a stream of variant candidates, and figures out where to break
        them up into windows and then gives that to the haplotype generator.

        Yields a dictionary.

        The basic procedure is as follows: loop through all variants. If the next variant is outside the
        window, then either extend the window, or, if we already have enough variants in the current
        window, then yield this window, and start a new one which will contain the next variant.

        The maximum allowed number of variants in a particular window is currently set to 3.
        Once we hit the end of the variants, we catch the StopIteration, and return whatever is in the the variant buffer.
        """
        windowVars = self.getWindowVariants(chromosome, start, end, sortedVariants, options.rlen, options.maxVariants)

        nVar = len(sortedVariants)
        nSnp = len([v for v in sortedVariants if abs(v.nAdded - v.nRemoved) == 0])

        logger.debug("There are %s vars in total in this region" %(nVar))
        logger.debug("There are %s SNPs and %s indels" %(nSnp, nVar-nSnp))
        logger.debug("There will be %s windows used in this region" %(len(windowVars)))

        nPos = len(windowVars)

        for index,vars in enumerate(windowVars):

            if options.callOnlyIndels:

                # Dirty dirty hack for 1kg indel calls
                foundIndel = False

                for v in vars:
                    if v.nAdded != v.nRemoved:
                        foundIndel = True
                        break

                if not foundIndel:
                    continue

            winStart = max(vars[0].minPos,start)
            winEnd = min(vars[-1].maxPos,end)

            if winStart == winEnd:
                if winEnd == end:
                    winStart = winStart - 1
                else:
                    winEnd = winEnd +1

            yield dict(chromosome=chromosome,startPos=winStart, endPos=winEnd, variants=vars)
