#cython: boundscheck=False
#cython: cdivision=True
#cython: nonecheck=False

"""
module containing various classes and functions for use in generating
and processing haplotypes.
"""

from __future__ import division

import cython
import logging
import math

cimport cython
cimport variant
cimport calign
cimport fastafile
cimport bamfileutils
cimport samtoolsWrapper
cimport cerrormodel

from calign cimport alignNoTraceback
from fastafile cimport FastaFile
from samtoolsWrapper cimport AlignedRead
from samtoolsWrapper cimport cAlignedRead
from variant cimport Variant

###################################################################################################

logger = logging.getLogger("Log")

###################################################################################################

cdef double PI = math.pi

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)
    double fabs(double)
    double sqrt(double)
    double pow(double, double)

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)
    void *memset(void *buffer, int ch, size_t count )

###################################################################################################

# Some nasty global variables
cdef list per_base_indel_errors = [2.9e-5, 2.9e-5, 2.9e-5, 2.9e-5, 4.3e-5, 1.1e-4, 2.4e-4, 5.7e-4, 1.0e-3, 1.4e-3] + [ 1.4e-3 + 4.3e-4*(n-10) for n in range(11,50) ]
cdef bytes homopolq = ''.join([chr(int(33.5 + 10*log( (idx+1)*q )/log(0.1) )) for idx,q in enumerate(per_base_indel_errors)])
cdef dict indel_error_model = {'nonrepetitive':'X',
                               1:'NKJHFA=854210/.-,,+**))(((\'\'\'&&&%%%$$$$#####"""""'}

###################################################################################################

cdef class Haplotype:
    """
    Class to encapsulate a single haplotype. This will store all the
    variants pertaining to this haplotype, as well as the reference
    sequence, all supporting reads, and the start and end positions of
    this haplotype in the reference.
    """
    def __init__(self, bytes refName, int startPos, int endPos, tuple variants, FastaFile refFile, int maxReadLength):
        """
        Constructor. Takes a tuple of variants and a
        fasta file of the referene sequence.
        Variants are sorted on instantiation
        """
        self.refName = refName
        self.refFile = refFile
        self.variants = variants
        self.varSeqLen = 0
        self.indelErrorModel = indel_error_model
        self.useIndelErrorModel = False

        cdef Variant v
        cdef int totalLengthChange = 0

        for v in variants:
            self.varSeqLen += (abs(v.nAdded - v.nRemoved))
            totalLengthChange += (v.nAdded - v.nRemoved)

        self.startPos = max(0, startPos - 2*maxReadLength)
        self.endPos = min(endPos + 2*maxReadLength, self.refFile.refs[self.refName].SeqLength)

        self.referenceSequence = self.refFile.getSequence(self.refName, self.startPos, self.endPos)
        self.hash = -1
        self.haplotypeSequence = None
        self.haplotypeSequence = self.getMutatedSequence()
        self.cHaplotypeSequence = self.haplotypeSequence
        self.hapLen = len(self.cHaplotypeSequence)

        if self.useIndelErrorModel:
            #self.localGapOpenQ = bytes("I"*len(self.haplotypeSequence))
            self.localGapOpenQ = cerrormodel.annotate_sequence(self.haplotypeSequence, self.indelErrorModel, 0, 24)
            self.cLocalGapOpenQ = self.localGapOpenQ

        self.cHomopolQ = homopolq
        self.likesByIndCache = NULL
        self.lenCache = 0

    def __dealloc__(self):
        """
        Clean up cache.
        """
        cdef int index = 0

        if self.lenCache > 0:

            for index from 0 <= index < self.lenCache:
                if self.likesByIndCache[index] != NULL:
                    free(self.likesByIndCache[index])

            free(self.likesByIndCache)

    def __copy__(self):
        """
        Make sure this never gets called for haplotypes.
        """
        raise StandardError, "Oh no! The bridge is gone!"

    def __cmp__(self, Haplotype other):
        """
        Comparison operator. This is implemented so that the haplotype
        instances may be sorted or stored in a priority queue.
        """
        return cmp(self.refName, other.refName) or cmp(self.startPos, other.startPos) or cmp(self.endPos, other.endPos)

    def __richcmp__(Haplotype self, Haplotype other, int opCode):
        """
        Comparison function:

        Are two haplotypes equal? Only return true if the mutated
        sequences are exactly equal.
        """
        # <
        if opCode == 0:
            if self.refName < other.refName:
                return True
            elif self.refName == other.refName and self.startPos < other.startPos:
                return True
            else:
                return False
        # <=
        elif opCode == 1:
            if self.refName > other.refName:
                return False
            elif self.refName == other.refName and self.startPos > other.startPos:
                return False
            else:
                return True
        # >
        elif opCode == 4:
            if self.refName > other.refName:
                return True
            elif self.refName == other.refName and self.startPos > other.startPos:
                return True
            else:
                return False
        # >=
        elif opCode == 5:
            if self.refName < other.refName:
                return False
            elif self.refName == other.refName and self.startPos < other.startPos:
                return False
            else:
                return True
        # ==
        if opCode == 2:
            if self.refName != other.refName:
                return False
            elif self.startPos != other.startPos:
                return False
            else:
                thisSeq = self.haplotypeSequence
                otherSeq = other.haplotypeSequence
                return thisSeq == otherSeq
        # !=
        elif opCode == 3:
            if self.refName != other.refName:
                return True
            elif self.startPos != other.startPos:
                return True
            else:
                thisSeq = self.haplotypeSequence
                otherSeq = other.haplotypeSequence
                return thisSeq != otherSeq
        else:
            raise StandardError, "Op code %s not implemented in haplotype__richcmp__()" %(opCode)

    def __hash__(self):
        """
        Implementing this function allows haplotypes to be hashed, and so stored in
        a set or dictionary. The supporting reads are not included in the hashing, as
        we want two haplotypes to give the same hash id if they have the same positions
        and sequences.
        """
        if self.hash == -1:
            self.hash = hash((self.refName, self.startPos, self.endPos, self.haplotypeSequence))

        return self.hash

    cdef int valid(self):
        """
        Check if this is a valid haplotype. If the variants overlap then we
        can't make a haplotype. Also optionally check if the standard form variant lies
        in the region provided - need to check at this level, rather than the variant level
        since this may depend on the reference sequence, and on the other variants in
        the haplotype
        """
        cdef int nVariants = len(self.variants)

        # Reference haplotype is always valid. And a single variant can't really overlap
        # with anything.
        if nVariants <= 1:
            return True

        cdef int index = 0
        cdef Variant thisVar
        cdef Variant nextVar

        # Check the normal variant starts and ends. The variants are sorted by co-ordinate, so if
        # the end of variant[i-1] occurs after the beginning of variant[i], then something is wrong.
        # Note, that the end of variant v is actually v.refPos + len(v.removed) - len(v.added).
        for index from 0 <= index < nVariants:

            thisVar = self.variants[index]

            # The end of the last variant can't overlap with anything.
            if index + 1 == nVariants:
                break

            nextVar = self.variants[index+1]

            # This should never happen
            if thisVar.minPos > nextVar.minPos:
                logger.error("Variants %s and %s are out of order. This should never happen." %(thisVar, nextVar))
                raise StandardError, "Variants out of order in haplotype!"

            if thisVar.maxPos >= nextVar.minPos:

                # Plain insertion can overlap, as this does not touch next reference base. But don't allow two
                # indels to overlap.
                if (thisVar.nAdded > 0 and thisVar.nRemoved == 0) and (nextVar.nAdded == 0 and nextVar.nRemoved == 0):
                    continue
                else:
                    #logger.debug("Overlapping Variants in haplotype %s. Rejecting." %(self))
                    return False

        # If we get here, then all is well.
        return True

    cdef double* alignReads(self, int individualIndex, cAlignedRead** start, cAlignedRead** end, int nReads, int nIndividuals):
        """
        """
        cdef int readIndex = 0
        cdef int cacheIndex = 0
        cdef double score = 0.0

        # Need to malloc cache array
        if self.lenCache == 0:

            self.likesByIndCache = <double**>(malloc(nIndividuals*sizeof(double*)))
            self.lenCache = nIndividuals

            for cacheIndex from 0 <= cacheIndex < self.lenCache:
                self.likesByIndCache[cacheIndex] = NULL

        # Not done.Calculate and cache values
        if self.likesByIndCache[individualIndex] == NULL:
            self.likesByIndCache[individualIndex] = <double*>(malloc(nReads*sizeof(double)))

            while start != end:

                if self.useIndelErrorModel:
                    score = alignNoTraceback(start[0], self.cHaplotypeSequence, self.startPos, self.varSeqLen, self.cLocalGapOpenQ, 3, 2, self.hapLen, False)
                else:
                    score = alignNoTraceback(start[0], self.cHaplotypeSequence, self.startPos, self.varSeqLen, self.cHomopolQ, 3, 2, self.hapLen, True)

                self.likesByIndCache[individualIndex][readIndex] = score
                start += 1
                readIndex += 1

        return self.likesByIndCache[individualIndex]

    cdef char* getReferenceSequence(self, prefix = 0):
        """
        Return the refernece sequence for the region covered by this haplotype. pretty shows where the variants are.
        """
        if prefix == 0 and self.referenceSequence != None:
            return self.referenceSequence

        seqMax = self.refFile.refs[self.refName].SeqLength
        self.referenceSequence = self.refFile.getSequence(self.refName, max(0, self.startPos - prefix), min(self.endPos + prefix, seqMax))
        return self.referenceSequence

    cdef char* getMutatedSequence(self):
        """
        Return the reference sequence mutated with all the variants being
        considered in this haplotype.

        Just to remind ourselves: SNPs are reported at the index of the changed
        base (0-indexed internally, and 1-indexed in VCF). Insertions are reported
        such that the index is that of the last reference base before the insertion.
        Deletions are reported such that the index is that of the first deleted base.
        """
        cdef Variant v
        cdef Variant firstVar

        if self.haplotypeSequence is None:

            if len(self.variants) == 0:
                self.haplotypeSequence = self.getReferenceSequence()
            else:
                currentPos = self.startPos

                # Get sequence up to one base before the first variant
                firstVar = self.variants[0]
                bitsOfMutatedSeq = [self.refFile.getSequence(self.refName, currentPos, firstVar.refPos)]
                currentPos = firstVar.refPos

                for v in self.variants:

                    # Move up to one base before the next variant, if we're not already there.
                    if v.refPos > currentPos:
                        bitsOfMutatedSeq.append(self.refFile.getSequence(self.refName, currentPos, v.refPos))
                        currentPos = v.refPos

                    nAdded = v.nAdded
                    nRemoved = v.nRemoved

                    # SNP/Mult-SNP
                    if nAdded == nRemoved:
                        bitsOfMutatedSeq.append(v.added)
                        currentPos += v.nAdded

                    # Deletion
                    elif nAdded == 0 and nRemoved > 0:

                        if v.refPos == currentPos:
                            bitsOfMutatedSeq.append(self.refFile.getCharacter(self.refName, v.refPos))
                            currentPos += 1

                        currentPos += nRemoved

                    # Insertion
                    elif nAdded > 0 and nRemoved == 0:

                        if v.refPos == currentPos:
                            bitsOfMutatedSeq.append(self.refFile.getCharacter(self.refName, v.refPos))
                            currentPos += 1

                        bitsOfMutatedSeq.append(v.added)

                    else:
                        logger.error("Found wonky variant in getMutatedSequence")
                        logger.error("Variants are %s" %(self.variants))
                        logger.error(v)
                        raise StandardError, "Variant type %s not yet supported in getMutatedSequence()" %(v)

                # Is this ok when currentPos == endPos?
                if currentPos < self.endPos:
                    bitsOfMutatedSeq.append(self.refFile.getSequence(self.refName, currentPos, self.endPos) )

                self.haplotypeSequence = ''.join(bitsOfMutatedSeq)

        return self.haplotypeSequence

    cdef list homopolymerLengths( self ):
        """
        Return a list of homopolymer lengths for the sequence surrounding
        each variant in this haplotype.
        """
        cdef tuple vsf = self.variants
        if len( vsf ) == 0:
            return []
        else:
            return [ self.homopolymerLengthForOneVariant( v ) for v in vsf ]

    cdef int homopolymerLengthForOneVariant(self, Variant variant):
        """
        Calculate and return the length of the largest homopolymer
        touching this variant. Compute homopolymer lengths on the
        left and right, and return the largest.
        """
        varChrom = variant.refName
        varPos = variant.refPos

        leftRefSeq = self.refFile.getSequence(varChrom, varPos-20, varPos)
        rightRefSeq = self.refFile.getSequence(varChrom, varPos+1, varPos + 21)

        if len(leftRefSeq) == 0 or len(rightRefSeq) == 0:
            return 0

        leftHpSize = 0
        rightHpSize = 0

        firstLeftChar = leftRefSeq[-1]
        firstRightChar = rightRefSeq[0]

        for char in reversed(leftRefSeq):
            if char == firstLeftChar:
                leftHpSize += 1
            else:
                break

        for char in rightRefSeq:
            if char == firstRightChar:
                rightHpSize += 1
            else:
                break

        if firstLeftChar != firstRightChar:
            return max(leftHpSize, rightHpSize)
        else:
            return leftHpSize + rightHpSize

    cdef bytes getSequenceContext(self, Variant variant):
        """
        Return the sequence surrounding this variant's position.
        """
        varChrom = variant.refName
        varPos = variant.refPos
        return self.refFile.getSequence(varChrom, varPos-10, varPos + 11)

    cdef dict vcfINFO(self):
        """
        Information to go in the vcf file INFO field - a two level dictionary of variant - field - value
        This can be augmented at the individual/population level with more information, or some of it can be
        removed before printing the output

        Include
        HP = Honopolymer tract length
        NR = Number of supporting reads
        CD = Coverage depth
        """
        cdef Variant variant
        cdef dict INFO = {}

        homopolymerLengths = self.homopolymerLengths()
        vsf = self.variants

        for varIndex, variant in enumerate( self.variants ):
            HP = homopolymerLengths[varIndex]
            TR = variant.nSupportingReads
            NF = variant.nFwdReads
            NR = variant.nRevReads
            SC = self.getSequenceContext(variant)

            INFO[vsf[varIndex]] = {'HP':[HP],'TR':[TR], 'NF':[NF], 'NR':[NR], 'SC':[SC]}

        return INFO

    cdef dict vcfFILTER(self):
        """
        Information to go in the FILTER field - a dictionary of variant to ( array of reasons to filter ).

        Includes
        hp10 - homopolymer length > 10
        ss - doesn't appear on reads covering both strands.
        """
        cdef dict FILTER = {}
        cdef Variant variant
        #homopolymerLengths = self.homopolymerLengths()
        cdef tuple vsf = self.variants

        for varIndex, variant in enumerate( self.variants ):
            varFILTER = []
        #    HP = homopolymerLengths[varIndex]

        #    if HP > 10:
        #        varFILTER.append('hp10')

            FILTER[vsf[varIndex]] = varFILTER
        return FILTER

    cdef double strandBiasPValue(self, Variant variant):
        """
        Calculate a binomial p-value for the strand bias
        """
        nReads = variant.nSupportingReads
        nFwdReads = variant.nFwdReads
        pVal = binomial(nFwdReads, nReads, 0.5)
        return pVal

    def __str__(self):
        """
        Generate a string representation of the haplotype. This is useful for
        debugging.
        """
        if len(self.variants) == 0:
            return '  Haplotype(*Reference*)  '

        vars = [str(v) for v in self.variants]
        string = "  Haplotype(" + ",".join(vars) + ")  "

        return string

    def __repr__(self):
        """
        The representation function. Called when printing the screen.
        """
        return self.__str__()

###################################################################################################

cdef double logFactorial(int x):
    """
    Return the logarithm of the factorial of x. Uses Stirling's approximation for large
    values.
    """
    cdef double ans = 0
    cdef int i = 0
    cdef double y = x

    if x < 15:
        for i from 1 <= i <= x:
            ans += log(i)
        return ans
    else:
        ans = y*log(y) + log(2.0*PI*y)/2 - y + (pow(y, -1))/12 - (pow(y,-3))/360 + (pow(y,-5))/1260 - (pow(y,-7))/1680 + (pow(y,-9))/1188
        return ans

###################################################################################################

cdef double binomial(int x, int size, double prob):
    """
    Optimised binomial probability function.
    """
    if x == size and prob == 1:
        return 1
    elif x != size and prob == 1:
        return 0
    elif x == 0 and prob == 0:
        return 1
    elif x == 0 and prob == 1:
        return 0

    cdef double logBinomCoefficient = logFactorial(size) - (logFactorial(x) + logFactorial(size-x))
    cdef double logBinomProb = x*log(prob) + (size-x)*log(1.0-prob)

    return exp(logBinomCoefficient + logBinomProb)

###################################################################################################

cdef Variant convertVariantToStandardFormat(Variant variant, FastaFile refFile, int maxReadLength):
    """
    Internally we use a representation which makes gives the reference position to be the base before the
    variant, and doesn't have a particular position in homopolymer runs. Here we convert to the standard
    format used in e.g. 1000g, For insertions we report the base after the variant, and report variants at the start
    of homopolymer runs. For deletions, we report the first deleted base, again starting at the beginning of homopolymer runs.

    To do this, we do a greedy match backwards on the mutated and reference sequences to work out where the
    variant will go.
    """
    cdef int nAdded = variant.nAdded
    cdef int nRemoved = variant.nRemoved

    # For SNPs, and multi-nucleotide substitutions, we don't need to do anything
    if nAdded == nRemoved:
        return variant

    # Hack. Leave variants alone if they occur at the extreme left end of the contigs, otherwise
    # we run out of sequence.
    if variant.refPos < 100:
        return variant

    # We look at the reference for this many bases either side of the variant, so that we can adjust its position
    # accordingly, if there are homopolymers etc.
    cdef int window = max(nAdded, nRemoved) + 10000
    cdef int seqMax = refFile.refs[variant.refName].SeqLength
    cdef int windowMin = max(1, variant.refPos - window)
    cdef int windowMax = min(variant.refPos + window, seqMax)

    cdef Haplotype hap = Haplotype(variant.refName, windowMin, windowMax, (variant,), refFile, 0)

    cdef bytes bytesRefSeq = hap.referenceSequence
    cdef bytes bytesHapSeq = hap.haplotypeSequence

    cdef char* refSequence = bytesRefSeq
    cdef char* hapSequence = hap.cHaplotypeSequence

    cdef int lenRef = len(refSequence)
    cdef int lenHap = hap.hapLen

    cdef int index = 0
    cdef int minLenRefHap = 0

    if lenRef < lenHap:
        minLenRefHap = lenRef
    else:
        minLenRefHap = lenHap

    # First loop forwards to find out how far to the right we can push it
    for index from 0 <= index < minLenRefHap:
        if hapSequence[index] != refSequence[index]:
            break;

    cdef int maxPos = variant.refPos - window + index + nAdded + nRemoved
    cdef int newPos = -1

    cdef char refChar
    cdef char hapChar

    cdef int hapIndex = 0
    cdef int refIndex = 0

    cdef bytes newAdded = bytes("")
    cdef bytes newRemoved = bytes("")

    cdef int effecctSize = 0
    cdef int insStart = 0
    cdef int delStart = 0
    cdef int lenNewAdded = 0
    cdef int lenNewRemoved = 0

    cdef Variant newVar

    # Loop backwards through the mutated and reference sequences
    for index from 0 <= index < minLenRefHap:

        hapIndex = (lenHap - index) - 1
        refIndex = (lenRef - index) - 1

        refChar = refSequence[refIndex]
        hapChar = hapSequence[hapIndex]

        # If the sequences are the same at this point, keep going.
        # Once we hit a difference, use that position.
        if hapChar != refChar:

            if nAdded > 0:
                newPos = windowMin + (lenRef - 1) - index
                insStart = newPos - windowMin + 1 # Position of first inserted base in haplotype sequence string
                newAdded = hapSequence[insStart: insStart + nAdded]

            if nRemoved > 0:
                newPos = windowMin + (lenRef - 1) - index - (nRemoved)
                delStart = newPos - windowMin + 1 # Position of first deleted base in reference sequence
                newRemoved = refSequence[delStart: delStart + nRemoved]

            if newPos > variant.refPos:
                logger.error("Old pos = %s new pos = %s" %(variant.refPos, newPos))
                logger.error(variant)
                logger.error(refSequence)
                logger.error(hapSequence)

            # Length of sequence on one side of this variant, which is affected by the variant. We set
            # this to be the difference in the max and min posiitions of the variant, + the length of inserted
            # or deleted sequence
            effectSize = len(newRemoved) + len(newAdded) + abs(newPos - variant.refPos) # For normalised variants, still consider reads covering original position

            newVar = Variant(variant.refName, newPos, newRemoved, newAdded, variant.nSupportingReads, variant.nFwdReads, 0, 0, 0)
            newVar.maxPos = newPos + effectSize
            newVar.minPos = newPos
            newVar.sumSqReadPos = variant.sumSqReadPos
            newVar.rmsReadPos = variant.rmsReadPos
            newVar.nFwdReadsWeighted = variant.nFwdReadsWeighted
            newVar.nRevReadsWeighted = variant.nRevReadsWeighted

            lenNewAdded = len(newAdded)
            lenNewRemoved = len(newRemoved)

            if lenNewAdded != nAdded or lenNewRemoved != nRemoved:
                logger.error("New variant in standard format %s is broken" %(newVar))
                logger.error("Original variant was %s" %(variant))
                logger.error(refSequence)
                logger.error(hapSequence)
                raise StandardError, "Error in variant conversion to standard format"

            return newVar

    # If we get to here, and we haven't returned, then it's almost certainly an error
    logger.error("Screw-up in converting variants to standard format")
    logger.error(refSequence)
    logger.error(hapSequence)
    logger.error(variant)

    raise StandardError, "Error in variant conversion to standard format. Variant = %s" %(variant)

###################################################################################################
