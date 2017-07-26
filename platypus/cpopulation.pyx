from __future__ import division

import logging
import platypus.vcf
import re
import math

from platypus.variant cimport Variant
from platypus.chaplotype cimport Haplotype
from platypus.cgenotype cimport DiploidGenotype
from platypus.samtoolsWrapper cimport AlignedRead
from platypus.samtoolsWrapper cimport cAlignedRead
from platypus.cwindow cimport bamReadBuffer

logger = logging.getLogger("Log")

ctypedef long long size_t
cdef double PI = math.pi
cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten

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


cdef extern from "string.h":
    char *strcpy(char *dest, char *src)
    long strlen(char *s)


cdef inline double double_max(double a, double b): return a if a >= b else b
cdef inline double double_min(double a, double b): return a if a <= b else b
cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b


vcfInfoSignature = {
    "FR":platypus.vcf.FORMAT('FR',1,0,'Float','Estimated population frequency',-1),
    "PP":platypus.vcf.FORMAT('PP',1,0,'Float','Posterior probability (phred scaled) that this variant segregates',-1),
    "TC":platypus.vcf.FORMAT('TC',1,0,'Integer','Total coverage at this locus',-1),
    "TCR":platypus.vcf.FORMAT('TCR',1,0,'Integer','Total reverse strand coverage at this locus',-1),
    "TCF":platypus.vcf.FORMAT('TCF',1,0,'Integer','Total forward strand coverage at this locus',-1),
    "TR":platypus.vcf.FORMAT('TR',1,0,'Integer','Total number of reads containing this variant',-1),
    "NF":platypus.vcf.FORMAT('NF',1,0,'Integer','Total number of forward reads containing this variant',-1),
    "NR":platypus.vcf.FORMAT('NR',1,0,'Integer','Total number of reverse reads containing this variant',-1),
    "SC":platypus.vcf.FORMAT('SC',1,1,'String','Genomic sequence 10 bases either side of variant position',-1),
    "RMP":platypus.vcf.FORMAT('RMP',1,0,'Float','RMS Position in reads of Variant',-1),
    "HP":platypus.vcf.FORMAT('HP',1,1,'Integer','Homopolmer run length',-1),
    "FPV":platypus.vcf.FORMAT('FPV',1,0,'Float','Forward strand p-value',-1),
    "RPV":platypus.vcf.FORMAT('RPV',1,0,'Float','Reverse strand p-value',-1),
    "ABPV":platypus.vcf.FORMAT('ABPV',1,0,'Float','Allele-bias p-value. Testing for low variant coverage',-1),
    "MMLQ":platypus.vcf.FORMAT('MMLQ',1,0,'Float','Median minimum base quality for bases around variant',-1),
}

vcfFilterSignature = {
    #"pp10":platypus.vcf.FORMAT('pp10',1,0,'Flag','Posterior probability phred-score is less than 10','.'),
    #"pp5":platypus.vcf.FORMAT('pp5',1,0,'Flag','Posterior probability phred-score is less than 5','.'),
    "ab":platypus.vcf.FORMAT('ab',1,0,'Flag','Variant fails allele-bias filter','.'),
    "sb":platypus.vcf.FORMAT('sb',1,0,'Flag','Variant fails strand-bias filter','.'),
    "badReads":platypus.vcf.FORMAT('badReads',1,0,'Flag','Variant supported only by reads with low quality bases close to variant position, and not present on both strands.','.'),
    "hp10":platypus.vcf.FORMAT('hp10',1,0,'Flag','Flanking sequence contains homopolymer of length 10 or greater','.'),
}


vcfFormatSignature = {
    "GT":platypus.vcf.FORMAT('GT',1,1,'String','Unphased genotypes','.'),
    "GL":platypus.vcf.FORMAT('GL',1,'.','Float','Genotype log-likelihoods (log10) for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites','.'),
    "GQ":platypus.vcf.FORMAT('GQ',1,1,'Integer','Genotype quality, as Phred score','.'),
    "NR":platypus.vcf.FORMAT('NR',1,1,'Integer','Number of reads covering variant in this sample','.'),
}


def getDataAndSampleNames(options):
    """
    Function to sort out the bam-file options and labels, for a dataset
    from a population.
    """
    logger.debug("Calling variants in population. Sorting out bam-files and labels")

    samples = None
    bamFiles = options.bamFiles
    nInd = len(bamFiles)

    if options.bamFiles == None:
        raise StandardError('For calling populations, must specify bamFiles')
    else:
        logger.debug("Data is in %s bam file(s)" %(len(bamFiles)))

    if options.labels == None:
        samples = [x.split("/")[-1] for x in options.bamFiles]
    else:
        samples = [re.search(options.labels, x).group(0) for x in options.bamFiles]

    return bamFiles,samples,nInd


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


cdef double binomial(int x, int size, double prob):
    """
    Optimised binomial probability function.
    """
    if x == size and prob == 1:
        return 1.0
    elif x != size and prob == 1:
        return 0.0
    elif x == 0 and prob == 0:
        return 1.0
    elif x == 0 and prob == 1:
        return 0.0
    elif x == 0 and size == 0:
        return 1.0

    cdef double logBinomCoefficient = logFactorial(size) - (logFactorial(x) + logFactorial(size-x))
    cdef double logBinomProb = x*log(prob) + (size-x)*log(1.0-prob)

    return exp(logBinomCoefficient + logBinomProb)


cdef class Population:
    """
    The Population class represents a population of genotypes. Each genotype
    represents the set of haplotypes assigned to a single individual in the genomic
    region under consideration.
    """
    def __init__(self, list genotypeCalls, list genotypeLikelihoods, list haplotypes, int nInd, int nGenotypes, list readBuffers, int verbosity, options):
        """
        Constructor. Takes the list of genotypes assigned to this population.
        Plus, a list of haplotypes and corresponding frequencies.
        """
        self.genotypeCalls = genotypeCalls
        self.genotypeLikelihoods = genotypeLikelihoods
        self.nGenotypes = nGenotypes
        self.haplotypes = haplotypes
        self.cutoff = 0.0 # Only report frequencies > this.
        self.nHaplotypes = len(haplotypes)
        self.nIndividuals = nInd
        self.readBuffers = readBuffers
        self.verbosity = verbosity
        self.badReadsWindow = options.badReadsWindow
        self.badReadsThreshold = options.badReadsThreshold
        self.sbThreshold = options.sbThreshold
        self.abThreshold = options.abThreshold
        logger.debug("Constructing population. nGenotypes = %s. nHaplotypes = %s. Lower frequency cut-off = %s." %(len(genotypeCalls), len(haplotypes), self.cutoff))

    def __dealloc__(self):
        """
        Make sure to clean up memory correctly.
        """
        cdef int index = 0
        cdef int genIndex = 0

        for index from 0 <= index < self.nIndividuals:
            free(self.genotypePosteriors[index])

        for genIndex from 0 <= genIndex < self.nGenotypes:
            free(self.haplotypeIndexes[genIndex])

        free(self.haplotypeIndexes)
        free(self.genotypePosteriors)
        free(self.frequencies)
        free(self.nReads)

    cdef dict vcfINFO(self):
        """
        INFO field for vcf output, two level dictionary of variant - info field - value
        """
        if self.vcfInfo is not None:
            return self.vcfInfo

        cdef Haplotype haplotype
        cdef DiploidGenotype genotype
        cdef Variant var
        cdef Variant variant
        cdef AlignedRead read
        cdef str PP
        cdef dict INFO = {}
        cdef dict hapInfo
        cdef double FR = 0.0
        cdef double posteriorProbThisVar = 0.0
        cdef double doubleFreq = 0.0
        cdef double weightedTCRVal = 0.0
        cdef double weightedTCFVal = 0.0
        cdef double fwdSbPVal = 0.0
        cdef double revSbPVal = 0.0
        cdef double abPVal = 0.0
        cdef int index = 0
        cdef int varPos = 0
        cdef int readStart = 0
        cdef int readEnd = 0
        cdef int readLen = 0
        cdef int lenVarAdded = 0
        cdef int TCVal = 0
        cdef int TCRVal = 0
        cdef int TCFVal = 0
        cdef int alleleCount = 0
        cdef int hapIndex = 0
        cdef int TN = 0
        cdef int minBaseQualInWindow = 0
        cdef int medMinBaseQual = 0
        cdef int windowStart = 0
        cdef int windowSize = self.badReadsWindow
        cdef int windowEnd = 0
        cdef int windowIndex = 0

        cdef list listOfMinBaseQuals = []
        cdef char* varAdded = NULL

        cdef cAlignedRead** pStartRead
        cdef cAlignedRead** pEndRead
        cdef cAlignedRead* pRead
        cdef bamReadBuffer readBuffer

        for hapIndex from 0 <= hapIndex < self.nHaplotypes:

            haplotype = self.haplotypes[hapIndex]
            hapINFO = haplotype.vcfINFO()

            for var,value in hapINFO.iteritems():

                if not INFO.has_key(var):
                    posteriorProbThisVar = self.calculatePosterior(var)
                    #logger.debug("Posterior for var %s = %s" %(var,posteriorProbThisVar))
                    PP = "%.0f" %(-10.0 * log10(double_max(1e-20, 1.0-posteriorProbThisVar)))
                    TR = var.nSupportingReads
                    NF = var.nFwdReads
                    NR = var.nRevReads

                    if NR + NF != TR:
                        logger.error("Broken variant %s. NR = %s. NF = %s. TR = %s" %(var, NR, NF, TR))

                    FR = self.frequencies[hapIndex]
                    RMP = var.rmsReadPos
                    INFO[var] = dict(HP=value['HP'], PP=[PP], TR=[TR], NR=[NR], NF=[NF],FR=[FR], SC=value['SC'], RMP=[round(RMP, 2)])
                else:
                    INFO[var]['FR'][0] += self.frequencies[hapIndex]

        for variant in INFO:

            TCVal = 0
            TCRVal = 0
            TCFVal = 0
            weightedTCRVal = 0.0
            weightedTFRVal = 0.0
            lenVarAdded = variant.nAdded
            varPos = variant.refPos
            varAdded = variant.added

            for index,genotype in enumerate(self.genotypeCalls):

                readBuffer = self.readBuffers[index]
                pStartRead = readBuffer.windowStart
                pEndRead = readBuffer.windowEnd

                while pStartRead != pEndRead:

                    pRead = pStartRead[0]
                    readLen = pRead.rlen
                    readStart = pRead.pos
                    readEnd = pRead.end

                    if readStart <= varPos and readEnd >= varPos:

                        #
                        if len(variant.added) != 1:
                            pStartRead += 1
                            continue

                        if pRead.seq[varPos - readStart] != varAdded[0]:
                            #print "Here. read = %s. var = %s" %(pRead.seq[varPos - readStart], varAdded[0])
                            pStartRead += 1
                            continue

                        # Loop over window around variant position, and store smallest quality
                        windowStart = max(0, varPos - readStart - (windowSize-1)//2)

                        # windowEnd is always 1 past the end
                        windowEnd = min(readLen, varPos - readStart + (windowSize-1)//2)

                        minBaseQualInWindow = 100

                        for windowIndex in range(windowStart, windowEnd):
                            minBaseQualInWindow = min(minBaseQualInWindow, pRead.qual[windowIndex])

                        listOfMinBaseQuals.append(minBaseQualInWindow)

                    pStartRead += 1

            listOfMinBaseQuals.sort()

            if len(listOfMinBaseQuals) > 0:
                medMinBaseQual = listOfMinBaseQuals[ len(listOfMinBaseQuals) // 2 ]
            else:
                medMinBaseQual = 100

            # Do the same loop over samples and reads again
            for index,genotype in enumerate(self.genotypeCalls):

                readBuffer = self.readBuffers[index]
                pStartRead = readBuffer.windowStart
                pEndRead = readBuffer.windowEnd

                while pStartRead != pEndRead:

                    pRead = pStartRead[0]
                    readLen = pRead.rlen
                    readStart = pRead.pos
                    readEnd = pRead.end

                    if readStart <= varPos and readEnd >= varPos:

                        TCVal += 1

                        if pRead.isReverse:
                            TCRVal += 1
                            weightedTCRVal += (1.0 - exp(mLTOT*(pRead.mapq)))
                        else:
                            TCFVal += 1
                            weightedTCFVal += (1.0 - exp(mLTOT*(pRead.mapq)))

                    pStartRead += 1


            if len(variant.added) == 1 and len(variant.removed) == 1:
                INFO[variant]['MMLQ'] = [medMinBaseQual]
            else:
                pass
                #INFO[variant]['MMLQ'] = []


            INFO[variant]['TC'] = [TCVal]
            INFO[variant]['TCR'] = [TCRVal]
            INFO[variant]['TCF'] = [TCFVal]
            doubleFreq = INFO[variant]['FR'][0]
            INFO[variant]['FR'][0] = "%1.4f" %(doubleFreq)

            fwdSbPVal = 0.0
            revSbPVal = 0.0
            abPVal = 0.0

            # I don't care if there are too many variant reads
            if ((variant.nFwdReadsWeighted + variant.nRevReadsWeighted) / <double>(TCVal)) >= (doubleFreq/2.0):
                abPVal = 1.0

            # Compare total variant reads to total coverage, and computed frequency
            else:
                for index in range(round(variant.nFwdReadsWeighted + variant.nRevReadsWeighted)+1):
                    abPVal += binomial(index, TCVal, doubleFreq/2.0)

            for index in range(round(variant.nFwdReadsWeighted)+1):
                fwdSbPVal += binomial(index, int(round(weightedTCFVal)), min(0.99,doubleFreq))

            for index in range(round(variant.nRevReadsWeighted)+1):
                revSbPVal += binomial(index, int(round(weightedTCRVal)), min(0.99,doubleFreq))

            if self.verbosity >= 3:
                logger.debug("Fwd pval = %s. rev pval = %s" %(fwdSbPVal, revSbPVal))

            INFO[variant]['FPV'] = ["%1.2e" %(fwdSbPVal)]
            INFO[variant]['RPV'] = ["%1.2e" %(revSbPVal)]
            INFO[variant]['ABPV'] = ["%1.2e" %(abPVal)]

        self.vcfInfo = INFO
        return INFO

    cdef dict vcfFILTER(self):
        """
        FILTER field for vcf output, two level dictionary of variant - info field - value
        """
        if self.vcfFilter is not None:
            return self.vcfFilter

        cdef DiploidGenotype genotype
        cdef Haplotype haplotype
        cdef Variant v
        cdef dict FILTER = {}
        cdef dict hapFILTER
        cdef dict INFO
        cdef dict infoThisVar
        cdef int nSamplesNoData = 0
        cdef int index = 0
        cdef int nVars = 0
        cdef int nVarsSB = 0
        cdef int nVarsAB = 0
        cdef int totalForwardVarReads = 0
        cdef int totalReverseVarReads = 0
        cdef int totalVarReads = 0
        cdef int totalReads = 0
        cdef int totalFwdReads = 0
        cdef int totalRevReads = 0
        cdef int varPP = 0
        cdef int nVarsFailingMMLQFilter = 0
        cdef int medMinQualBases = 0
        cdef double varFreq = 0.0
        cdef double fwdPVal = 0.0
        cdef double revPVal = 0.0
        cdef double abPVal = 0.0

        for index,genotype in enumerate(self.genotypeCalls):
            if genotype is None:
                nSamplesNoData += 1

        for haplotype in self.haplotypes:
            hapFILTER = haplotype.vcfFILTER()
            FILTER.update(hapFILTER)

        INFO = self.vcfINFO()

        for v in FILTER:

            infoThisVar = INFO[v]

            varFreq = float(infoThisVar['FR'][0])
            totalForwardVarReads = int(infoThisVar['NF'][0])
            totalReverseVarReads = int(infoThisVar['NR'][0])
            totalVarReads = int(infoThisVar['TR'][0])
            totalReads = int(infoThisVar['TC'][0])
            totalFwdReads = int(infoThisVar['TCF'][0])
            totalRevReads = int(infoThisVar['TCR'][0])
            fwdPVal = float(infoThisVar['FPV'][0])
            revPVal = float(infoThisVar['RPV'][0])
            abPVal = float(infoThisVar['ABPV'][0])
            varPP = int(infoThisVar['PP'][0])
            medMinQualBases = int(infoThisVar.get('MMLQ', [100])[0])

            if medMinQualBases < self.badReadsThreshold and not (totalForwardVarReads > 0 and totalReverseVarReads > 0):
                nVarsFailingMMLQFilter += 1

            nVars += 1

            if totalForwardVarReads == 0 and totalFwdReads != 0:
                if fwdPVal < self.sbThreshold:
                    nVarsSB += 1

            if totalReverseVarReads == 0 and totalRevReads != 0:
                if revPVal < self.sbThreshold:
                    nVarsSB += 1

            if abPVal < self.abThreshold:
                nVarsAB += 1

        if nVarsFailingMMLQFilter == nVars:
            FILTER[v].append('badReads')

        if nVarsSB == nVars:
            FILTER[v].append('sb')

        if nVarsAB == nVars:
            FILTER[v].append('ab')

        self.vcfFilter = FILTER
        return FILTER

    cdef dict vcfVARIANTS(self):
        """
        dict of position - array of variants
        """
        if self.vcfVariants is not None:
            return self.vcfVariants

        cdef Haplotype haplotype
        cdef DiploidGenotype gen
        cdef Variant v
        cdef dict variantDict = {}
        cdef int i = 0

        for i,haplotype in enumerate(self.haplotypes):

            if self.frequencies[i] < self.cutoff:
                continue

            hapCalled = False

            for gen in self.genotypeCalls:
                if gen and haplotype in gen.haplotypes:
                    hapCalled = True

            if not hapCalled:
                continue

            for v in haplotype.variants:
                try:
                    variantDict[v.refPos].add(v)
                except:
                    variantDict[v.refPos] = set([v])

        self.vcfVariants = variantDict
        return variantDict

    cdef double calculatePosterior(self, Variant var):
        """
        Calculate for ecah variant the posterior probability that that variant (in standard format)
        is segregating in the population
        """
        cdef Haplotype hap
        cdef tuple vsf
        cdef double* genotypeLikelihoodsThisIndividual
        cdef double freqsPrimeByHapIndex[1000]
        cdef double prior = var.calculatePrior()
        cdef double totalProb = 0.0
        cdef double sumFreqs = 0.0
        cdef double thisFreq = 0.0
        cdef double ratio = 0.0
        cdef double likelihoodThisGenotype = 0.0
        cdef double variantPosterior = 0.0
        cdef double sumProbVariantThisIndividual = 0.0
        cdef double sumProbNoVariantThisIndividual = 0.0
        cdef double sumLogProbVariant = 0.0
        cdef double sumLogProbNoVariant = 0.0
        cdef double r_prime = 0.0
        cdef double s_prime = 0.0
        cdef double factor = 0.0
        cdef int i = 0
        cdef int j = 0
        cdef int r = 0
        cdef int s = 0
        cdef int nReadsThisInd = 0
        cdef int logOfMinFloat = -708

        # TODO: This restriction should be enforced earlier.
        if self.nHaplotypes > 1000:
            raise StandardError, "Too many haplotypes (%s). This should never happen" %(self.nHaplotypes)

        if self.verbosity >= 3:
            logger.debug("N haplotypes = %s. n Ind = %s" %(self.nHaplotypes, self.nIndividuals))

        # Re-scale haplotype frequencies for those not containing this variant, and
        # store them.
        for i from 0 <= i < self.nHaplotypes:

            if self.verbosity >= 3:
                logger.debug("Freq for hap %s = %s" %(self.haplotypes[i], self.frequencies[i]))

            hap = self.haplotypes[i]
            vsf = hap.variants

            if var not in vsf:
                freqsPrimeByHapIndex[i] = self.frequencies[i]
                sumFreqs += self.frequencies[i]
            else:
                freqsPrimeByHapIndex[i] = 0.0

        if self.verbosity >= 3:
            logger.debug("Sum freqs = %s" %(sumFreqs))

        if sumFreqs > 0:
            for i from 0 <= i < self.nHaplotypes:
                if freqsPrimeByHapIndex[i] == 0.0:
                    continue
                else:

                    if self.verbosity >= 3:
                        logger.debug("before scaling, freqsPrimeByHapIndex[%s] = %s" %(i, freqsPrimeByHapIndex[i]))

                    freqsPrimeByHapIndex[i] /= sumFreqs

                    if self.verbosity >= 3:
                        logger.debug("after scaling, freqsPrimeByHapIndex[%s] = %s" %(i, freqsPrimeByHapIndex[i]))

        # Compute the variant posterior.
        for i from 0 <= i < self.nIndividuals:

            nReadsThisInd = self.nReads[i]

            if nReadsThisInd == 0:
                continue

            genotypeLikelihoodsThisIndividual = self.genotypePosteriors[i]
            sumProbVariantThisIndividual = 0.0
            sumProbNoVariantThisIndividual = 0.0

            for j from 0 <= j < self.nGenotypes:

                likelihoodThisGenotype = genotypeLikelihoodsThisIndividual[j]

                r = self.haplotypeIndexes[j][0]
                s = self.haplotypeIndexes[j][1]

                factor = 1.0

                if r != s:
                    factor = 2.0

                sumProbVariantThisIndividual += (factor * self.frequencies[r] * self.frequencies[s] * likelihoodThisGenotype)

                # r_prime and s_prime will be 0.0 for all haplotypes containing the variant
                r_prime = freqsPrimeByHapIndex[r]
                s_prime = freqsPrimeByHapIndex[s]
                sumProbNoVariantThisIndividual += (factor * r_prime * s_prime * likelihoodThisGenotype)

            if sumProbVariantThisIndividual > 0:
                sumLogProbVariant += log(sumProbVariantThisIndividual)
            else:
                sumLogProbVariant += logOfMinFloat

            if sumProbNoVariantThisIndividual > 0:
                sumLogProbNoVariant += log(sumProbNoVariantThisIndividual)
            else:
                sumLogProbNoVariant += logOfMinFloat

        ratio = exp(sumLogProbNoVariant - sumLogProbVariant)
        variantPosterior = prior / (prior + ratio * (1.0 - prior))

        if self.verbosity >= 3:
            logger.debug("For variant %s, Prior = %s. Ratio = %s. Post = %s. SumLogProbNoVar = %s. SumLogProbVar = %s" %(var,prior,ratio,variantPosterior,sumLogProbNoVariant,sumLogProbVariant))

        return variantPosterior


cdef class Caller:
    """
    The Population Caller takes a list of haplotypes and a list of reads and
    computes a MLE for the population frequency of each of the variants.
    """
    def __init__(self, options):
        """
        Constructor does nothing.
        """
        self.options = options

    def __dealloc__(self):
        """
        De-allocate memory. Make sure to do this correctly for 2-d
        arrays!
        """
        cdef int index = 0

        for index from 0 <= index < self.nIndividuals:
            free(self.cisr[index])

        free(self.cisr)

    cdef setup(self, list haplotypes, list genotypes, int nInd, int ploidy, int verbosity, list readBuffers):
        """
        Constructor for the population caller. The initial loop over reads, calculating
        likelihoods for each read, given a particular haplotype, is done here.

        We cache posterior values, the indexes in 'haplotypes' of haplotypes for each genotype,
        which are used in each step of the EM algorithm but do not change.
        """
        assert ploidy == 2

        self.ploidy = ploidy
        self.nGenotypes = len(genotypes)
        self.haplotypes = haplotypes
        self.genotypes = genotypes
        self.nHaplotypes = len(haplotypes)
        self.nIndividuals = nInd
        self.readBuffers = readBuffers
        self.verbosity = verbosity

        self.genotypePosteriors = <double**>(calloc(self.nIndividuals, sizeof(double*)))
        self.cisr = <double**>(calloc(self.nIndividuals, sizeof(double*)))
        self.haplotypeIndexes = <int**>(malloc(self.nGenotypes*sizeof(int*)))
        self.nReads = <int*>(malloc(self.nIndividuals*sizeof(int)))

        cdef int genotypeIndex = 0
        cdef int individualIndex = 0
        cdef int haplotypeIndex
        cdef int index = 0
        cdef double logLikelihood = 0.0
        cdef double likelihood = 0.0

        cdef dict hapIndDict = dict([(hap,index) for index,hap in enumerate(self.haplotypes)])
        cdef DiploidGenotype genotype
        cdef Haplotype hap1
        cdef Haplotype hap2
        cdef int nReadsThisInd = 0
        cdef bamReadBuffer theBuffer

        cdef double maxLogLike = -1e7

        for individualIndex from 0 <= individualIndex < self.nIndividuals:

            self.genotypePosteriors[individualIndex] = <double*>(calloc(self.nGenotypes, sizeof(double)))
            self.cisr[individualIndex] = <double*>(calloc(self.nGenotypes, sizeof(double)))
            theBuffer = self.readBuffers[individualIndex]
            nReadsThisInd = theBuffer.windowEnd - theBuffer.windowStart
            self.nReads[individualIndex] = nReadsThisInd

            for genotypeIndex from 0 <= genotypeIndex < self.nGenotypes:

                genotype = genotypes[genotypeIndex]

                if individualIndex == 0:
                    self.haplotypeIndexes[genotypeIndex] = <int*>(calloc(ploidy, sizeof(int)))
                    hap1 = genotype.hap1
                    hap2 = genotype.hap2
                    self.haplotypeIndexes[genotypeIndex][0] = <int>(hapIndDict[hap1])
                    self.haplotypeIndexes[genotypeIndex][1] = <int>(hapIndDict[hap2])

                if nReadsThisInd == 0:
                    self.genotypePosteriors[individualIndex][genotypeIndex] = 0.0
                else:
                    logLikelihood = genotype.calculateDataLikelihood(theBuffer.windowStart, theBuffer.windowEnd, individualIndex, nReadsThisInd, self.nIndividuals)

                    if logLikelihood > maxLogLike:
                        maxLogLike = logLikelihood

                    self.genotypePosteriors[individualIndex][genotypeIndex] = logLikelihood

                    # DEBUG OUTPUT
                    if verbosity >= 3:
                        logger.debug("LL for genotype %s is %s. L = %s. Number of reads = %s" %(genotype, logLikelihood, exp(logLikelihood), self.nReads[individualIndex]))

        # Normalising.
        for individualIndex from 0 <= individualIndex < self.nIndividuals:
            for genotypeIndex from 0 <= genotypeIndex < self.nGenotypes:
                self.genotypePosteriors[individualIndex][genotypeIndex] = exp(self.genotypePosteriors[individualIndex][genotypeIndex] - maxLogLike)

                # DEBUG OUTPUT
                if verbosity >= 3:
                    genotype = genotypes[genotypeIndex]
                    logger.debug("Normalised L for genotype %s is %s. Max log like = %s" %(genotype, self.genotypePosteriors[individualIndex][genotypeIndex], maxLogLike))

    cdef Population call(self, int maxIters = 100):
        """
        Estimate the haplotype frequencies and genotypes for this population and locus.
        Terminate the algorithm once the maximum change across all components of the
        frequency is < eps, or after maxIters iterations.
        """
        cdef list genotypeCalls = []
        cdef list genotypeLikelihoods = []
        cdef double* csr
        cdef double* freq = <double*>(calloc(self.nHaplotypes, sizeof(double)))
        cdef double eps = min(1e-3, 1.0 / (self.nIndividuals*self.ploidy*2))
        cdef double maxChange = eps + 1
        cdef double uniformFreq = 1.0/self.nHaplotypes
        cdef double maxCSR = -1.0
        cdef double thisCSR = 0.0
        cdef int iters = 0
        cdef int index = 0
        cdef int genotypeIndex = 0
        cdef int indexOfMaxCSR = -1
        cdef int nReadsThisInd = 0

        # Uniform initial frequency
        for index from 0 <= index < self.nHaplotypes:
            freq[index] = uniformFreq

        while maxChange > eps and iters < maxIters:
            maxChange = self.EMiteration(freq)
            iters += 1

        logger.debug("MaxChange = %s. nIterations = %s" %(maxChange, iters))

        # Now choose the most likely genotypes
        for index from 0 <= index < self.nIndividuals:

            nReadsThisInd = self.nReads[index]

            if nReadsThisInd == 0:
                genotypeCalls.append(None)
                genotypeLikelihoods.append(0.0)
                continue

            csr = self.cisr[index]
            indexOfMaxCSR = -1
            maxCSR = 0.0

            # Find genotype for this individual with best csr value
            for genotypeIndex from 0 <= genotypeIndex < self.nGenotypes:

                thisCSR = csr[genotypeIndex]

                if indexOfMaxCSR == -1 or thisCSR > maxCSR:
                    maxCSR = thisCSR
                    indexOfMaxCSR = genotypeIndex

            gt = self.genotypes[indexOfMaxCSR]
            genotypeCalls.append(gt)
            genotypeLikelihoods.append(maxCSR)

        pop = Population(genotypeCalls, genotypeLikelihoods, self.haplotypes, self.nIndividuals, self.nGenotypes, self.readBuffers, self.verbosity, self.options)
        pop.frequencies = freq
        pop.genotypePosteriors = self.genotypePosteriors
        pop.haplotypeIndexes = self.haplotypeIndexes
        pop.nReads = self.nReads

        return pop

    cdef double EMiteration(self, double* freq):
        """
        Perform one EM update, given initial frequency freq. Returns the updated frequency
        """
        cdef double* genotypePosteriorsThisIndividual
        cdef double* csr
        cdef double maxChange = 0.0
        cdef double csrSum = 0.0
        cdef double tmp = 0.0
        cdef double thisCSR = 0.0
        cdef double oldFreq = 0.0
        cdef double freqChange = 0.0
        cdef double absFreqChange = 0.0
        cdef double newFreq = 0.0
        cdef double posteriorThisGenotype = 0.0
        cdef int i = 0
        cdef int j = 0
        cdef int k = 0
        cdef int s = 0
        cdef int r = 0
        cdef int nIndWithData = 0
        cdef int nReadsThisInd = 0

        for i from 0 <= i < self.nIndividuals:

            nReadsThisInd = self.nReads[i]

            if nReadsThisInd == 0:
                continue

            genotypePosteriorsThisIndividual = self.genotypePosteriors[i]
            csr = self.cisr[i]
            csrSum = 0.0
            nIndWithData += 1

            for j from 0 <= j < self.nGenotypes:

                posteriorThisGenotype = genotypePosteriorsThisIndividual[j]
                s = self.haplotypeIndexes[j][0]
                r = self.haplotypeIndexes[j][1]
                thisCSR = posteriorThisGenotype*freq[s]*freq[r]*(1 + (r != s))
                csr[j] = thisCSR
                csrSum += thisCSR

            # Normalisation
            if csrSum > 0.0:
                for j from 0 <= j < self.nGenotypes:
                    csr[j] /= csrSum

        for k from 0 <= k < self.nHaplotypes:

            tmp = 0.0

            for i from 0 <= i < self.nIndividuals:

                nReadsThisInd = self.nReads[i]

                if nReadsThisInd == 0:
                    continue

                csr = self.cisr[i]

                for j from 0 <= j < self.nGenotypes:

                    thisCSR = csr[j]
                    s = self.haplotypeIndexes[j][0]
                    r = self.haplotypeIndexes[j][1]

                    if k == s or k == r:

                        if s == r:
                            tmp += (2.0*thisCSR)
                        else:
                            tmp += (thisCSR)

            oldFreq = freq[k]
            newFreq = tmp/(self.ploidy*nIndWithData)
            freq[k] = newFreq
            freqChange = fabs(oldFreq - newFreq)

            if freqChange > maxChange:
                maxChange = freqChange

        return maxChange
