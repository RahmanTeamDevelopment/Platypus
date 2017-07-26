import logging

from platypus.samtoolsWrapper cimport AlignedRead
from platypus.samtoolsWrapper cimport cAlignedRead
from platypus.fastafile cimport FastaFile

logger = logging.getLogger("Log")

cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)
    double fabs(double)
    double sqrt(double)


cdef class Variant(object):
    """
    Class to encapsulate information for all common variant types. The basic
    idea is to tread all variants as replacements, i.e. a certain sequence of
    bases is removed from the reference, and a certain sequence is added. A SNP is
    a replacement of 1 base with 1 different base, a deletion is a replacement of one
    sequence with a smaller sequence (or none at all), and an insertion is a replacement
    of one sequence with a longer sequence.
    """
    def __init__(self, bytes refName, int refPos, char* removed, char* added, int nSupportingReads, int nFwdReads, int mapQ, int baseQ, int readPos):
        """
        Constructor. Takes the name of the reference sequence (e.g. 'chr1'),
        the position of the variant (given by the left-most co-ordinate of the
        changed sequence (so a SNP co-ordinate is the exact location of the SNP,
        an indel co-ordinate is the location of the first inserted base, and a deletion
        co-ordinate is the location of the first deleted base)), the removed and added
        sequences, and the supporting read.
        """
        # Make sure we don't go negative here.
        refPos = max(0, refPos)

        self.refName = refName
        self.minPos = refPos

        self.nAdded = len(added)
        self.nRemoved = len(removed)
        self.freqPrior = 0.0

        # SNP
        if self.nRemoved == 1 and self.nAdded == 1:
            self.maxPos = refPos #+ 1 # Want to always consider pairs of SNPs together, if they occur together.

        # Multi-nucleotide substitution
        elif self.nRemoved == self.nAdded:
            self.maxPos = refPos + (self.nAdded - 1)

        # Indel/complex variation
        else:
            self.maxPos = refPos + self.nRemoved + self.nAdded

        self.refPos = refPos
        self.removed = removed
        self.added = added
        self.nSupportingReads = nSupportingReads
        self.nFwdReads = nFwdReads
        self.nRevReads = self.nSupportingReads - self.nFwdReads
        self.hashValue = -1
        self.sumSqReadPos = (readPos*readPos)

        if self.nSupportingReads > 0:
            self.rmsReadPos = sqrt(self.sumSqReadPos/self.nSupportingReads)

        if nFwdReads == 1:
            self.nFwdReadsWeighted = (1.0 - exp(mLTOT*baseQ))
            self.nRevReadsWeighted = 0.0
        else:
            self.nFwdReadsWeighted = 0.0
            self.nRevReadsWeighted = (1.0 - exp(mLTOT*baseQ))

    cdef double calculatePrior(self):
        """
        Calculate and return the prior probability for this
        variant.
        """
        if self.freqPrior != 0.0:
            return self.freqPrior

        cdef double prior = 0.0

        # Basic prior for SNPs is 1e-3, and for indels is 1e-4
        if self.nAdded == 1 and self.nRemoved == 1:
            prior = 1e-3

        # Multi-nucleotide substitution
        elif self.nAdded == self.nRemoved:
            prior = 5e-5*(0.1**self.nAdded)
        # Insertion
        elif self.nAdded > 0 and self.nRemoved == 0:
            #prior = 1e-4*(0.25**nAdded)
            prior = 1e-4*(0.33**self.nAdded)
        # Deletion
        elif self.nAdded == 0 and self.nRemoved > 0:
            prior = 1e-4*(0.8**self.nRemoved)
        else:
            raise StandardError, "Do not know how to calculate prior for variant %s." %(self.__str__())

        return prior

    cdef void addVariant(self, Variant other):
        """
        Add supporting data from another variant instance.
        """
        self.nSupportingReads += other.nSupportingReads
        self.nFwdReads += other.nFwdReads
        self.nRevReads += other.nRevReads
        self.nFwdReadsWeighted += other.nFwdReadsWeighted
        self.nRevReadsWeighted += other.nRevReadsWeighted
        self.sumSqReadPos += other.sumSqReadPos

        if self.nSupportingReads > 0:
            self.rmsReadPos = sqrt(self.sumSqReadPos/self.nSupportingReads)

    def __hash__(self):
        """
        Implementing this function allows variants to be hashed, and so stored in
        a set or dictionary. The supporting reads are not included in the hashing, as
        we want two variants to give the same hash id if they have the same position
        and added/removed sequences.
        """
        if self.hashValue == -1:
            self.hashValue = hash( (self.refName, self.refPos, self.removed, self.added) )

        return self.hashValue

    def __cmp__(self, Variant other):
        """
        Comparison operator. This is implemented so that the variant
        instances may be sorted or stored in a priority queue. Variants
        are compared first by position, and then by numbers of added and removed
        bases. For variants at the same position, SNPs should always come first, the
        insertions/deletions.
        """
        return cmp(self.refName, other.refName) or cmp(self.refPos, other.refPos) or cmp(abs(self.nAdded - self.nRemoved), abs(other.nAdded - other.nRemoved))

    def __richcmp__(Variant self, Variant other, int opCode):
        """
        Comparison function:

        Are two variants equal? Only return true if the positions and added and
        removed sequences are exactly equal.
        """
        cdef int thisRefPos = self.refPos
        cdef int otherRefPos = other.refPos
        cdef int thisNAdded = self.nAdded
        cdef int otherNAdded = other.nAdded
        cdef int thisNRemoved = self.nRemoved
        cdef int otherNRemoved = other.nRemoved
        cdef int thisLenChange = abs(thisNAdded - thisNRemoved)
        cdef int otherLenChange = abs(otherNAdded - otherNRemoved)
        cdef bytes thisRefName = self.refName
        cdef bytes otherRefName = other.refName
        cdef bytes thisAdded = self.added
        cdef bytes otherAdded = other.added
        cdef bytes thisRemoved = self.removed
        cdef bytes otherRemoved = other.removed

        # <
        if opCode == 0:
            if thisRefName < otherRefName:
                return True
            elif thisRefName == otherRefName and thisRefPos < otherRefPos:
                return True
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisLenChange < otherLenChange:
                return True
            else:
                return False
        # <=
        elif opCode == 1:
            if thisRefName > otherRefName:
                return False
            elif thisRefName == otherRefName and thisRefPos > otherRefPos:
                return False
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisLenChange > otherLenChange:
                return False
            else:
                return True
        # ==
        elif opCode == 2:
            if thisRefName != otherRefName:
                return False
            elif thisRefPos != otherRefPos:
                return False
            elif thisRemoved != otherRemoved:
                return False
            elif thisAdded != otherAdded:
                return False
            else:
                return True
        # !=
        elif opCode == 3:
            if thisRefName != otherRefName:
                return True
            elif thisRefPos != otherRefPos:
                return True
            elif thisRemoved != otherRemoved:
                return True
            elif thisAdded != otherAdded:
                return True
            else:
                return False
        # >
        elif opCode == 4:
            if thisRefName > otherRefName:
                return True
            elif thisRefName == otherRefName and thisRefPos > otherRefPos:
                return True
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisLenChange > otherLenChange:
                return True
            else:
                return False
        # >=
        elif opCode == 5:
            if thisRefName < otherRefName:
                return False
            elif thisRefName == otherRefName and thisRefPos < otherRefPos:
                return False
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisLenChange < otherLenChange:
                return False
            else:
                return True

    def __str__(self):
        """
        Another way of printing the variant.
        """
        string = "Var(%s:%s" %(self.refName,self.refPos)

        if self.nRemoved > 0:
            string += (' -' + self.removed)
        if self.nAdded > 0:
            string += (' +' + self.added)

        string += " nReads = %s, nFwd = %s, nRev = %s" %(self.nSupportingReads,self.nFwdReads,self.nRevReads)
        string += ")"
        return string

    def __repr__(self):
        """
        __repr__ is called when you do "print variant", and will return a short string
        describing the variant.
        """
        return str(self)


cdef class VariantCandidateGenerator(object):
    """
    A class to generate variant candidates from a bunch of reads.
    """
    def __init__(self, tuple region, list bamFiles, FastaFile referenceFile, int minMapQual, int minFlank, int minBaseQual, int maxCoverage, int maxReadLength, int verbosity=2, int genSNPs=1, int genIndels=1):
        """
        Constructor. Takes read buffer, reader, reference, and a set of quality arguments for the variant
        candidates. Create a storage splace for variant candidates, and store the values of some flags which
        are used in the pysam CIGAR information (these should really be module level variables).
        """
        self.CIGAR_M = 0 # Match
        self.CIGAR_I = 1 # Insertion
        self.CIGAR_D = 2 # Deletion
        self.CIGAR_N = 3 # Skipped region from reference
        self.CIGAR_S = 4 # Soft clipping. Sequence is present in read
        self.CIGAR_H = 5 # Hard clipping. Sequence is not present in read
        self.CIGAR_P = 6 # Padding. Used for padded alignment

        self.minMapQual = minMapQual
        self.minBaseQual = minBaseQual
        self.minFlank = minFlank
        self.variantHeap = {} # List of variants
        self.refFile = referenceFile
        self.bamFiles = bamFiles
        self.rname = region[0]
        self.refSeqStart = max(1, region[1]-2000) # Don't try to fetch anything < 0
        self.refSeqEnd = min(region[2]+2000, self.refFile.refs[region[0]].SeqLength) # Don't try to fetch anything > seqLength
        self.pyRefSeq = self.refFile.getSequence(self.rname, self.refSeqStart, self.refSeqEnd) # Cache this
        self.refSeq = self.pyRefSeq
        self.rStart = region[1]
        self.rEnd = region[2]
        self.maxCoverage = maxCoverage
        self.maxReadLength = maxReadLength
        self.verbosity = verbosity
        self.genSNPs = genSNPs
        self.genIndels = genIndels

    cdef addVariantToList(self, Variant var):
        """
        Check if the variant is already on the heap; if so, increment the number of
        supporting reads; if not, add it.
        """
        cdef Variant theVar = self.variantHeap.get(var)

        if theVar is None:
            self.variantHeap[var] = var

            if self.verbosity >= 3:
                logger.debug("Adding new variant %s to candidate list" %(var))

        else:
            theVar.addVariant(var)

            if self.verbosity >= 3:
                logger.debug("Adding new variant %s to existing variant in candidate list" %(var))

    cdef getSnpCandidatesFromReadSegment(self, cAlignedRead* read, char* readSeq, char* readQual, int readStart, int readOffset, int refOffset, int lenSeqToCheck):
        """
        Get all SNP candidates from a particular read segement.

        Args:
        char* readSeq -- The complete sequence of bases in the read.
        char* readQual -- The complete sequence of quality scores in the read.
        int readStart -- The starting position, in the reference sequence, of the read.
        int readOffset -- If we're not starting from the beginning of the read, then how far along the read sequence to start.
        int refOffset -- If the read contains indels, we need to offset our reference position accordingly, by this much.
        int lenSeqToCheck -- The number of bases to check.
        """
        cdef int index = 0
        cdef int baseQual = 0
        cdef int readLength = read.rlen
        cdef int mapQ = read.mapq
        cdef int isReverse = read.isReverse
        cdef int nContiguousSNPs = 0
        cdef int readIndex = 0
        cdef int refIndex = 0
        cdef char readChar
        cdef char refChar

        cdef int mSNPRefStart = 0
        cdef int mSNPReadStart = 0
        cdef bytes mSNPRef
        cdef bytes mSNPAlt

        for index from 0 <= index < lenSeqToCheck:

            readIndex = index + readOffset
            refIndex = (index + refOffset + readStart) - self.refSeqStart
            readChar = readSeq[readIndex]
            refChar = self.refSeq[refIndex]
            baseQual = readQual[readIndex]

            if readChar != refChar and readChar != 'N' and refChar != 'N' and baseQual >= self.minBaseQual:
                nContiguousSNPs += 1

                if isReverse:
                    self.addVariantToList(Variant(self.rname, refIndex + self.refSeqStart, <bytes>refChar, <bytes>readChar, 1, 0, mapQ, baseQual, readIndex+1))
                else:
                    self.addVariantToList(Variant(self.rname, refIndex + self.refSeqStart, <bytes>refChar, <bytes>readChar, 1, 1, mapQ, baseQual, readIndex+1))

                # Possible multi-nucleotide substitution
                #if nContiguousSNPs > 1:
                #    mSNPRefStart = refIndex - (nContiguousSNPs- 1)
                #    mSNPReadStart = readIndex - (nContiguousSNPs - 1)
                #    mSNPRef = self.refSeq[mSNPRefStart: mSNPRefStart+nContiguousSNPs]
                #    mSNPAlt = readSeq[mSNPReadStart: mSNPReadStart+nContiguousSNPs]

                #    if isReverse:
                #        self.addVariantToList(Variant(self.rname, mSNPRefStart + self.refSeqStart, mSNPRef, mSNPAlt, 1, 0), 0)
                #    else:
                #        self.addVariantToList(Variant(self.rname, mSNPRefStart + self.refSeqStart, mSNPRef, mSNPAlt, 1, 1), 1)
            else:
                nContiguousSNPs = 0

    cdef getVariantCandidatesFromSingleRead(self, cAlignedRead* read):
        """
        Check a single read for variant candidates. Variants are flagged by the CIGAR string. Pysam reports
        the CIGAR string information as a list of tuples, where each tuple is a pair, and the first element
        gives the type of feature (match, insertion or deletion), and the second element gives the number of
        nucleotides associated. For example, [(0, 1), (1, 2), (0, 1)] is a 1 base match, a 2 base insertion,
        and a 1 base match.
        """
        cdef int readStart = read.pos
        cdef int readLength = read.rlen
        cdef int mapQ = read.mapq
        cdef int flag = 0
        cdef int length = 0
        cdef int refOffset = 0
        cdef int readOffset = 0
        cdef int cigarIndex = 0
        cdef int cigarLength = read.cigarLen
        cdef int isReverse = read.isReverse
        cdef char* readQual = read.qual
        cdef char* readSeq = read.seq
        cdef bytes insertedSequence = None
        cdef bytes deletedSequence = None

        for cigarIndex from 0 <= cigarIndex < cigarLength:

            flag = read.cigarOps[cigarIndex]
            length = read.cigarLens[cigarIndex]

            # An insertion take us further along the read, but not the reference
            if flag == self.CIGAR_I and self.genIndels:

                # Skip this insertion if it isn't flanked by a matching sequence >= minFlank
                if cigarIndex > 0 and read.cigarOps[cigarIndex-1] == self.CIGAR_M and read.cigarLens[cigarIndex-1] >= self.minFlank:
                    pass
                elif cigarIndex < cigarLength-1 and read.cigarOps[cigarIndex+1] == self.CIGAR_M and read.cigarLens[cigarIndex+1] >= self.minFlank:
                    pass
                else:
                    readOffset += length
                    continue

                insertedSequence = readSeq[readOffset : readOffset+length]

                if insertedSequence.count("N") == 0 and self.genIndels:
                    if isReverse:
                        self.addVariantToList(Variant(self.rname, readStart+refOffset-1, "", insertedSequence, 1, 0, mapQ, readQual[readOffset], readOffset+1))
                    else:
                        self.addVariantToList(Variant(self.rname, readStart+refOffset-1, "", insertedSequence, 1, 1, mapQ, readQual[readOffset], readOffset+1))

                readOffset += length

            # A deletion take us further along the reference, but not the read
            elif flag == self.CIGAR_D:

                # Skip this deletion if it isn't flanked by a matching sequence >= minFlank
                if cigarIndex > 0 and read.cigarOps[cigarIndex-1] == self.CIGAR_M and read.cigarLens[cigarIndex-1] >= self.minFlank:
                    pass
                elif cigarIndex < cigarLength-1 and read.cigarOps[cigarIndex+1] == self.CIGAR_M and read.cigarLens[cigarIndex+1] >= self.minFlank:
                    pass
                else:
                    refOffset += length
                    continue

                deletedSequence = self.refFile.getSequence(self.rname, readStart+refOffset, (readStart+refOffset+length))

                # Don't look at deletions with Ns in them
                if deletedSequence.count("N") == 0 and self.genIndels:
                    if isReverse:
                        self.addVariantToList(Variant(self.rname, readStart+refOffset-1, deletedSequence, "", 1, 0, mapQ, readQual[readOffset], readOffset+1))
                    else:
                        self.addVariantToList(Variant(self.rname, readStart+refOffset-1, deletedSequence, "", 1, 1, mapQ, readQual[readOffset], readOffset+1))

                refOffset += length

            # A match take us further along the reference and the read
            elif flag == self.CIGAR_M:

                # Don't generate SNP candidates from matching sequences < minFlank
                if length < self.minFlank:
                    readOffset += length
                    refOffset += length
                    continue

                if self.genSNPs:
                    self.getSnpCandidatesFromReadSegment(read, readSeq, readQual, readStart, readOffset, refOffset, length)

                readOffset += length
                refOffset += length

            # Skipped region from the reference.
            elif flag == self.CIGAR_N:
                readOffset += length
                refOffset += length

            # Soft clipping. Sequence is present in read, but we should ignore it.
            elif flag == self.CIGAR_S:
                readOffset += length

            # Hard clipping. Sequence is not present in read.
            elif flag == self.CIGAR_H:
                continue

            # Padding. We do nothing here.
            elif flag == self.CIGAR_P:
                continue

            # Other kinds of flag.
            else:
                continue

    cpdef getListOfSnpPositions(self, AlignedRead read):
        """
        Return a list of positions at which the read has SNP candidates.
        """
        cdef int readStart = read.pos()
        cdef int flag = 0
        cdef int length = 0
        cdef int refOffset = 0
        cdef int readOffset = 0
        cdef int cigarIndex = 0
        cdef int cigarLength = read.getCigarLength()
        cdef int index = 0
        cdef int baseQual = 0
        cdef int readIndex = 0
        cdef int refIndex = 0
        cdef int readLength = read.rlen()
        cdef int isReverse = read.is_reverse()
        cdef char* readQual = read.qual()
        cdef char* readSeq = read.seq()
        cdef char readChar
        cdef char refChar
        cdef list snpPos = []

        for cigarIndex from 0 <= cigarIndex < cigarLength:

            flag = read.getCigarOpCode(cigarIndex)
            length = read.getCigarOpLength(cigarIndex)

            # An insertion take us further along the read, but not the reference
            if flag == self.CIGAR_I:
                readOffset += length

            # A deletion take us further along the reference, but not the read
            elif flag == self.CIGAR_D:
                refOffset += length

            # A match take us further along the reference and the read
            elif flag == self.CIGAR_M:

                for index from 0 <= index < length:

                    readIndex = index + readOffset
                    refIndex = (index + refOffset + readStart) - self.refSeqStart
                    readChar = readSeq[readIndex]
                    refChar = self.refSeq[refIndex]
                    baseQual = readQual[readIndex]

                    if readChar != refChar and readChar != 'N' and refChar != 'N' and baseQual >= self.minBaseQual:
                        if isReverse:
                            snpPos.append([readLength - readIndex - 1])
                        else:
                            snpPos.append([readIndex])

                readOffset += length
                refOffset += length

            # Other kinds of CIGAR flag. Most likely padding information (MAQ uses this).
            else:
                continue

        if isReverse:
            snpPos.reverse()
        return snpPos

    cpdef getListOfIndelPositions(self, AlignedRead read):
        """
        Return a list of positions at which the read has indel candidates.
        """
        cdef int readStart = read.pos()
        cdef int flag = 0
        cdef int length = 0
        cdef int refOffset = 0
        cdef int readOffset = 0
        cdef int cigarIndex = 0
        cdef int cigarLength = read.getCigarLength()
        cdef int index = 0
        cdef int baseQual = 0
        cdef int readIndex = 0
        cdef int refIndex = 0
        cdef int readLength = read.rlen()
        cdef int isReverse = read.is_reverse()
        cdef char* readQual = read.qual()
        cdef char* readSeq = read.seq()
        cdef char readChar
        cdef char refChar
        cdef list indelPos = []

        for cigarIndex from 0 <= cigarIndex < cigarLength:

            flag = read.getCigarOpCode(cigarIndex)
            length = read.getCigarOpLength(cigarIndex)

            # An insertion take us further along the read, but not the reference
            if flag == self.CIGAR_I:

                if isReverse:
                    indelPos.append([readLength - (readOffset + length)])
                else:
                    indelPos.append([readOffset])

                readOffset += length

            # A deletion take us further along the reference, but not the read
            elif flag == self.CIGAR_D:

                if isReverse:
                    indelPos.append([readLength - (readOffset + length)])
                else:
                    indelPos.append([readOffset])

                refOffset += length


            # A match take us further along the reference and the read
            elif flag == self.CIGAR_M:
                readOffset += length
                refOffset += length

            # Other kinds of CIGAR flag. Most likely padding information (MAQ uses this).
            else:
                continue

        return indelPos

    cdef addCandidatesFromReads(self, cAlignedRead** readStart, cAlignedRead** readEnd):
        """
        Loop through all reads, and flag candidate variants.
        """
        cdef int nReads = 0

        while readStart != readEnd:
            self.getVariantCandidatesFromSingleRead(readStart[0])
            nReads += 1
            readStart += 1

            if nReads % 100000 == 0:
                logger.info("Detected %s variants in %s reads" %(len(self.variantHeap), nReads))

    cdef list getCandidates(self):
        return sorted(self.variantHeap.values())
