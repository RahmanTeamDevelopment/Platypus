from __future__ import division

import logging

from platypus.samtoolsWrapper cimport Samfile
from platypus.samtoolsWrapper cimport IteratorRow
from platypus.samtoolsWrapper cimport createRead
from platypus.samtoolsWrapper cimport destroyRead
from platypus.fastafile cimport FastaFile
from platypus.chaplotype cimport Haplotype
from platypus.cfilter cimport getFilteredHaplotypes


logger = logging.getLogger("Log")


cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)


cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)
    double fabs(double)
    int abs(int)


cdef int bisectReadsLeft(cAlignedRead** reads, int testPos, int nReads):
    """
    Specialisation of bisection algorithm for array of
    read pointers.
    """
    cdef int low = 0
    cdef int high = nReads
    cdef int mid = 0

    while low < high:

        mid = (low + high) / 2

        if reads[mid].pos < testPos:
            low = mid + 1
        else:
            high = mid

    return low


cdef int bisectReadsRight(cAlignedRead** reads, int testPos, int nReads):
    """
    Specialisation of bisection algorithm for array of
    read pointers.
    """
    cdef int low = 0
    cdef int high = nReads
    cdef int mid = 0

    while low < high:

        mid = (low + high) / 2

        if testPos < reads[mid].pos:
            high = mid
        else:
            low = mid + 1

    return low


cdef int checkAndTrimRead(cAlignedRead* theRead, cAlignedRead* theLastRead, int maxBadQualBases, int minGoodQualBases):
    """
    Performs various quality checks on the read, and trims read (i.e. set q-scores to zero). Returns
    true if read is ok, and false otherwise.
    """
    # Filter out low quality reads, i.e. those with < minGoodQualBases with Q scores >= 20
    if theRead.rlen - theRead.nQualsBelow20 < minGoodQualBases:
        return False

    # Filter out reads with too many bad bases
    if theRead.nQualsBelow20 > maxBadQualBases:
        return False

    # Remove unmapped reads
    if theRead.isUnMapped:
        return False

    # Remove duplicate reads
    if theRead.isDuplicate:
        return False

    # If the insert size is < read length, then we almost certainly have adapter contamination, in which case we should
    # trim both reads. For now, though, just skip these reads.
    #
    # Actually, skipping is probably the right thing to do, as some of these reads will be mapped to the wrong location
    # in the genome.
    if theRead.insertSize > 0 and (not theRead.isReverse) and theRead.mateIsReverse and theRead.insertSize < theRead.rlen:
        return False

    if theRead.insertSize < 0 and theRead.isReverse and (not theRead.mateIsReverse) and -1*theRead.insertSize < theRead.rlen:
        return False

    # Check if this read is actually a duplicate. TODO: store library tag and check.
    if theLastRead != NULL:
        if theRead.pos == theLastRead.pos and theRead.rlen == theLastRead.rlen:

            # For paired reads, check mate's position
            if theRead.isPaired:
                if theLastRead.matePos == theRead.matePos:
                    return False

            # For single reads, just check pos and length of reads
            else:
                return False

    ## Any read that gets passed here will be used, but low quality tails will be trimmed, and
    ## any overlapping pairs, from small fragments, will have the overlapping bits trimmed.
    cdef int index = 0

    # Trim low-quality tails
    for index from 1 <= index <= theRead.rlen:

        if theRead.qual[theRead.rlen - index] < 5:
            theRead.qual[theRead.rlen - index] = 0
        else:
            break

    # Trim low-quality heads
    for index from 0 <= index < theRead.rlen:

        if theRead.qual[index] < 5:
            theRead.qual[index] = 0
        else:
            break

    # Trim overlapping part of forward read, in pairs where the read length is greater than the insert size
    # N.B Insert size is from start of forward read to end of reverse read, i.e. fragment size. This is done to
    # remove duplicate information, which gives systematic errors when pcr errors have occured in library prep.

    if theRead.insertSize > 0 and (not theRead.isReverse) and theRead.mateIsReverse and theRead.insertSize < 2*theRead.rlen:

        for index from 1 <= index < (2*theRead.rlen - theRead.insertSize) + 1:
            theRead.qual[theRead.rlen - index] = 0

    return True


cdef class bamReadBuffer(object):
    """
    Utility class for bufffering reads from a single BAM file, so we only make a single pass
    through the data in each BAM in the loop through windows.
    """
    def __init__(self, str fileName, char* chrom, int start, int end, int minMapQual, int maxReads, int minGoodQualBases, int maxBadQualBases):
        """
        Constructor.
        """
        cdef int lenReadArray = 100

        self.chrom = chrom
        self.startBase = start
        self.endBase = end
        self.longestRead = 0
        self.cReads = <cAlignedRead**>malloc(lenReadArray*sizeof(cAlignedRead))

        cdef Samfile reader = Samfile(fileName, 'rb')
        cdef IteratorRow readIterator = reader.fetch(self.chrom, self.startBase, self.endBase)

        cdef int readStart = 0
        cdef int readEnd = 0
        cdef int readLength = 0
        cdef int currentIndel = 0
        cdef int readCount = 0
        cdef int mapQualRejectCount = 0
        cdef int baseQualRejectCount = 0
        cdef int readOk = 0

        cdef cAlignedRead** tempReadArray = NULL

        while readIterator.cnext():

            # Remove reads with low mapping qualities
            if readIterator.b.core.qual < minMapQual:
                mapQualRejectCount += 1
                continue

            if readCount >= lenReadArray:
                tempReadArray = <cAlignedRead**>(realloc(self.cReads, 2*sizeof(cAlignedRead*)*lenReadArray))

                if tempReadArray == NULL:
                    raise StandardError, "Could not re-allocate read array in read buffer"
                else:
                    self.cReads = tempReadArray
                    lenReadArray = 2*lenReadArray

            self.cReads[readCount] = createRead(readIterator.b)

            if readCount > 0:
                readOk = checkAndTrimRead(self.cReads[readCount], self.cReads[readCount - 1], maxBadQualBases, minGoodQualBases)
            else:
                readOk = checkAndTrimRead(self.cReads[readCount], NULL, maxBadQualBases, minGoodQualBases)

            if not readOk:
                destroyRead(self.cReads[readCount])
                continue

            if self.cReads[readCount].isUnMapped:
                logger.warning("Un-mapped read in window. Chrom = %s. Pos = %s." %(chrom,self.cReads[readCount].pos))

            readStart = self.cReads[readCount].pos
            readEnd = self.cReads[readCount].end
            readLength = readEnd - readStart

            if readLength > self.longestRead:
                self.longestRead = readLength

            readCount += 1

            if readCount > maxReads:
                logger.info("Too many reads (>%s) in buffer. Stopping now.")
                break

        self.nReads = readCount
        self.start = self.cReads
        self.end = self.cReads + self.nReads # One past the end
        logger.debug("Read %s reads from file %s." %(readCount,fileName))
        reader.close()

    def __dealloc__(self):
        """
        Clean up memory
        """
        cdef int index = 0

        for index from 0 <= index < self.nReads:
            destroyRead(self.cReads[index])

        free(self.cReads)

    cdef setWindowPointers(self, int start, int end):
        """
        Set the windowStart and windowEnd pointers to point to the first
        and last+1 reads covering this window.
        """
        if self.nReads == 0:
            self.windowStart = self.cReads
            self.windowEnd = self.cReads
            return

        cdef int firstOverlapStart = max(1, start - self.longestRead)
        cdef int startPosOfReads = bisectReadsLeft(self.cReads, firstOverlapStart, self.nReads)
        cdef int endPosOfReads = bisectReadsLeft(self.cReads, end, self.nReads)

        while startPosOfReads < self.nReads and self.cReads[startPosOfReads].end <= start:
            startPosOfReads += 1

        # Cython has no (*p) dereference operator. Usef p[0] instead.
        self.windowStart = self.cReads + startPosOfReads
        self.windowEnd = min(self.cReads + endPosOfReads, self.end)

        if startPosOfReads > endPosOfReads:
            logger.info("Start pos = %s. End pos = %s. Read start pos = %s. end pos = %s" %(start, end, startPosOfReads, endPosOfReads))
            logger.info("There are %s reads here." %(self.nReads))
            logger.info("Buffer start pos = %s. end pos = %s." %(self.startBase, self.endBase))
            raise StandardError, "This should never happen. Read start pointer > read end pointer!!"


cdef list getHaplotypesInWindow(dict window, int nReads, FastaFile refFile, int ploidy, int maxCoverage, int minMapQual, int minBaseQual, int maxHaplotypes, int maxVariants, int maxReadLength):
    """
    Takes a window with some potential variants. Using all data in this region, calculate what the best possible
    haplotyes are and return them. If there are zero reads covering this window, then simply return the reference
    haplotype; if the mean coverage is > maxCoverage, then raise an exception.

    Returns:

    haplotypes: a list of haplotypes
    """
    cdef str windowChr = window['chromosome']
    cdef int windowStart = window['startPos']
    cdef int windowEnd = window['endPos']
    cdef list variants = window['variants']
    cdef int nVar = len(variants)
    cdef int readLength = 0
    cdef Haplotype refHaplotype = Haplotype(windowChr, windowStart, windowEnd, (), refFile, maxReadLength)

    # Make sure that we log regions of zero coverage. This should not really happen.
    if nReads == 0:
        logger.warning("Coverage is zero in window %s:%s-%s. Variants considered in this window are %s" %(windowChr, windowStart, windowEnd, variants))
        return [refHaplotype]

    # We throw if coverage is too high, to avoid heavily duplicated
    # regions.
    coverage = (nReads * maxReadLength) / ( windowEnd-windowStart+maxReadLength)

    if coverage > maxCoverage:
        logger.warning("Coverage %s exceeds the upper limit (%s). There are %s reads." %(coverage, maxCoverage, nReads))
        return [refHaplotype]

    # Take the top 'maxHaplotypes' haplotypes, as ranked by likelihood
    return getFilteredHaplotypes(refHaplotype, variants, nVar, refFile, windowChr, windowStart, windowEnd, maxHaplotypes, maxReadLength)
