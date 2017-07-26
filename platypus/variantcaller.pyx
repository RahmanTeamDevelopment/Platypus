from __future__ import division

import window
import logging
import variantutils
import datetime
import time

import platypus.vcf
import platypus.cpopulation
cimport platypus.fastafile

from operator import attrgetter

from platypus.chaplotype cimport Haplotype
from platypus.cgenotype cimport generateAllGenotypesFromHaplotypeList
from platypus.cgenotype cimport DiploidGenotype
from platypus.cwindow cimport bamReadBuffer
from platypus.samtoolsWrapper cimport AlignedRead
from platypus.cpopulation cimport Population
from platypus.cpopulation cimport Caller
from platypus.variant cimport Variant
from platypus.variant cimport VariantCandidateGenerator
from platypus.fastafile cimport FastaFile
from platypus.variantFilter cimport filterVariants

logger = logging.getLogger("Log")
vcfHeader = [('fileDate',datetime.date.fromtimestamp(time.time())), ('source','Platypus_Version_0.1.5')]
nSupportingReadsGetter = attrgetter("nSupportingReads")


cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)


cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)


cdef callVariantsInSmallRegion(chrom, int start, int end, list bamFiles, FastaFile refFile, options, windowGenerator, outputFile, vcf, list samples):
    """
    """
    cdef bamReadBuffer theReadBuffer
    cdef AlignedRead read
    cdef Haplotype refHaplotype
    cdef Haplotype hap
    cdef VariantCandidateGenerator varCandGen
    cdef Variant v

    cdef list readBuffers = []
    cdef list allReadsByIndividual
    cdef list readsThisSample
    cdef list varHaps
    cdef list varGens
    cdef list varsByPost
    cdef list variants
    cdef list allVarHaplotypes
    cdef list allUniqueHaplotypes
    cdef list reads
    cdef list rawBamVariants
    cdef list filteredBamVariants

    cdef set variantSet = set()
    cdef dict thisWindow = None

    cdef int nReadsThisWindow = 0
    cdef int regionSize = 0
    cdef int nVar = 0
    cdef int index = 0

    cdef Population call
    cdef Caller caller

    # Buffer BAM data
    cdef int nReadsInBuffer = 0

    for fileName in bamFiles:
        theReadBuffer = bamReadBuffer(fileName, chrom, start, end, options.minMapQual,options.maxReads,options.minGoodQualBases,options.maxBadQualBases)
        nReadsInBuffer += theReadBuffer.nReads

        if nReadsInBuffer > options.maxReads:
            logger.warning("Too many reads (>%s) in region %s:%s-%s. Skipping this region." %(options.maxReads, chrom, start, end))
            return

        readBuffers.append(theReadBuffer)

    # Generate initial variant candidate list
    if options.getVariantsFromBAMs:
        varCandGen = VariantCandidateGenerator((chrom,start,end), bamFiles, refFile, options.minMapQual, options.minFlank, options.minBaseQual, options.maxReads, options.rlen, options.verbosity, options.genSNPs, options.genIndels)

        for theReadBuffer in readBuffers:
            varCandGen.addCandidatesFromReads(theReadBuffer.start, theReadBuffer.end)

        rawBamVariants = varCandGen.getCandidates()
        filteredBamVariants = filterVariants(rawBamVariants, refFile, options.rlen, options.minReads, options.strandFilter, options.maxSize, options.verbosity)
        variantSet.update(filteredBamVariants)

    if options.sourceFile:
        variantSource = variantutils.VariantCandidateReader(options.sourceFile)
        vcfFileVariants = [v for v in variantSource.Variants(chrom,start,end, options.freqAsPrior)]
        variantSet.update(vcfFileVariants)

    cdef list allSortedVariants = sorted(variantSet)

    if options.printVarsAndExit:
        for v in allSortedVariants:
            print v
        return

    # Generate windows and loop through
    for thisWindow in windowGenerator.WindowsAndVariants(chrom, start, end, allSortedVariants, options):

        variants = thisWindow["variants"]
        windowStart = thisWindow['startPos']
        windowEnd = thisWindow['endPos']
        windowGenerator.nWindows += 1

        if options.verbosity >= 3:
            for v in variants:
                logger.debug(v)

        if windowGenerator.nWindows % 250 == 0:
            logger.info("There are %s variants in window %s:%s-%s" %(len(variants), chrom, windowStart, windowEnd))

        refHaplotype = Haplotype(chrom, thisWindow['startPos'], thisWindow['endPos'], (), refFile, options.rlen)
        regionSize = windowEnd - windowStart
        nReadsThisWindow = 0
        nVar = len(variants)

        try:

            for index,theReadBuffer in enumerate(readBuffers):
                theReadBuffer.setWindowPointers(windowStart, windowEnd)
                nReadsThisWindow += (theReadBuffer.windowEnd - theReadBuffer.windowStart)

            # Filter variants by posterior
            if nVar > options.maxVariants:

                # For pathological cases, first sort the variants by coverage,
                # and take the top 50 variants.
                if nVar > 50:
                    logger.info("Taking top 50 variant by coverage")
                    variants = sorted(sorted(variants, key=nSupportingReadsGetter, reverse=True)[0:50])
                    nVar = len(variants)

                logger.info("There are %s variants. Filtering variants with EM method" %(nVar))
                varHaps = [refHaplotype] + [Haplotype(chrom, windowStart, windowEnd, (v,), refFile, options.rlen) for v in variants]
                varGens = generateAllGenotypesFromHaplotypeList(options.ploidy, varHaps)
                caller = cpopulation.Caller(options)
                caller.setup(varHaps, varGens, options.nInd, options.ploidy, options.verbosity, readBuffers)
                call = caller.call()
                varsByPost = []

                for v in variants:
                    post = call.calculatePosterior(v)
                    varsByPost.append( (post,v) )

                varsByPost.sort(reverse=True)

                variants = sorted([theTuple[1] for theTuple in varsByPost[:options.maxVariants]])
                thisWindow['variants'] = variants

            # Create haplotype list using all data from all samples. Always consider the reference haplotype.
            allVarHaplotypes = cwindow.getHaplotypesInWindow(thisWindow, nReadsThisWindow, refFile, options.ploidy, options.maxReads, options.minMapQual, options.minBaseQual, options.maxHaplotypes, options.maxVariants, options.rlen)
            allUniqueHaplotypes = list(set([refHaplotype] + allVarHaplotypes))
            nUniqueHaplotypes = len(allUniqueHaplotypes)

            if options.verbosity >= 3:
                for h in allUniqueHaplotypes:
                    logger.debug(h)

            # Create genotype list
            allGenotypes = generateAllGenotypesFromHaplotypeList(options.ploidy, allUniqueHaplotypes)

            if nUniqueHaplotypes <= 1:
                logger.debug("No haplotypes to check in window %s" %(thisWindow))
                continue

            caller = cpopulation.Caller(options)
            caller.setup(allUniqueHaplotypes, allGenotypes, options.nInd, options.ploidy, options.verbosity, readBuffers)
            call = caller.call()
            outputCallToVCF(call, vcf, samples, samples, refFile, outputFile, options)

        except Exception, e:
            logger.exception('Exception in window %d-%d.' %(thisWindow['startPos'],thisWindow['endPos']))
            continue

###################################################################################################

def generateGenotypesAndVariantsInAllRegionsAndHandleExceptions(region, fileName, options, lock):
    """
    Call and return variants and genotypes in all regions, ensuring that,
    when multiple processes are run, each region gets processed only once.
    """
    assert options.ploidy == 2
    assert options.bufferSize > 0
    assert options.bufferSize <= options.processRegionSize

    logger.debug("Constructing VariantCaller class with options %s" %(options))

    options.bamFiles, samples, options.nInd = platypus.cpopulation.getDataAndSampleNames(options)
    options = options

    chromosome,start,end = region

    # Cache reference sequence for this region. This should result in a substantial speed-up for
    # the fastafile.getSequence function.
    cdef FastaFile refFile = platypus.fastafile.FastaFile(options.refFile, options.refFile + ".fai")
    refFile.setCacheSequence(chromosome, start-25000, end+25000)

    outputFile = open(fileName, 'wb')

    vcf = platypus.vcf.VCF()
    vcf.setheader(vcfHeader + [('platypusOptions', str(options))])
    vcf.setsamples(samples)
    vcf.setinfo(platypus.cpopulation.vcfInfoSignature)
    vcf.setfilter(platypus.cpopulation.vcfFilterSignature)
    vcf.setformat(platypus.cpopulation.vcfFormatSignature)
    vcf.writeheader(outputFile)

    allBamFiles = options.bamFiles
    windowGenerator = window.WindowGenerator(refFile, options.sourceFile)
    windowGenerator.nWindows = 0

    for i in range(start, end, options.bufferSize):
        thisStart = i
        thisEnd = min(i+options.bufferSize, end)
        callVariantsInSmallRegion(chromosome, thisStart, thisEnd, allBamFiles, refFile, options, windowGenerator, outputFile, vcf, samples)

    outputFile.close()


def generateGenotypesAndVariantsInAllRegions(tuple args):
    """
    A wrapper for the above funcion, to catch exceptions.
    """
    (region, fileName, options, lock) = args

    try:
        generateGenotypesAndVariantsInAllRegionsAndHandleExceptions(region, fileName, options, lock)
    except Exception, e:
        logger.exception(e.message)
        raise e


cdef outputCallToVCF(Population call, vcfFile, list allSamples, list samplesThisPop, FastaFile refFile, outputFile, options):
    """
    Output a call to a vcf file. Uses the signature of the callingModule.
    the file is opened and the header written when this object is initialised.
    The call object must have vcfINFO,vcfFILTER, genotypes, and variants methods.
    variants is a dictionary of POS-[array of variants]. Genotypes is a list of
    genotypes in the same order as the samples header of the vcf file.
    this code prints for each of the variants, the appropriate line.
    """
    cdef DiploidGenotype AAGenotype
    cdef DiploidGenotype ABGenotype
    cdef DiploidGenotype BBGenotype
    cdef DiploidGenotype genotype

    cdef Haplotype testHap
    cdef Haplotype refHap
    cdef Haplotype varHap
    cdef Haplotype hapOne
    cdef Haplotype hapTwo

    cdef Variant theVar
    cdef bamReadBuffer theBuffer

    cdef tuple varsHapOne
    cdef tuple varsHapTwo

    cdef dict varDict = call.vcfVARIANTS()
    cdef dict info = call.vcfINFO()
    cdef dict filter = call.vcfFILTER()
    cdef dict lineinfo

    cdef list POSs = sorted(varDict.keys())
    cdef list readsThisIndividual
    cdef list genotypeLikelihoodsThisInd
    cdef list FR
    cdef list PP
    cdef list NF
    cdef list NR
    cdef list TR
    cdef list GT
    cdef list linefilter
    cdef list genotypeCalls = call.genotypeCalls
    cdef list genotypeLikelihoods = call.genotypeLikelihoods

    cdef int i = 0
    cdef int variantIndex = 0
    cdef int genotypeQuality = 0
    cdef int nIndividuals = call.nIndividuals
    cdef int nVariants = 0
    cdef int tnReadsThisIndividual = 0

    for POS in POSs:

        variants = sorted([v for v in list(varDict[POS]) if float(info[v]['PP'][0]) >= options.minPosterior])
        nVariants = len(variants)

        if options.verbosity >= 3:
            logger.debug("Haps are %s" %(call.haplotypes))
            logger.debug("Logging mutated sequences...")
            logger.debug("Lengths are %s" %([len(str(refHap.getMutatedSequence())) for refHap in call.haplotypes]))

            for refHap in call.haplotypes:
                logger.debug(str(refHap.getMutatedSequence())[50:-50])

            logger.debug("Unfiltered variants are %s" %(sorted([v for v in list(varDict[POS])])))
            logger.debug("Unfiltered variant info is %s" %([info[v] for v in sorted(list(varDict[POS]))]))
            logger.debug("Filtered Variants are %s" %(variants))
            logger.debug("NVariants = %s" %(nVariants))

        if nVariants == 0:
            continue

        chrom = variants[0].refName
        ref,alt = refAndAlt(chrom, POS, variants, refFile)
        id = "."            # For now - if you want to look it up, go ahead.
        qual = -1           # Needs to be added
        format = ['GT:GL:GQ:NR']
        linefilter = []
        lineinfo = info[variants[0]]
        FR = []             # Can be generalised to all fields of variable length
        PP = []
        NF = []
        NR = []
        TR = []


        for var in variants:

            #logger.debug("Info for variant %s is %s" %(var, info[var]))
            linefilter.extend( [f for f in filter[var] if ( f in vcfFile.getfilter() )] )

            freq = info[var]['FR']
            posterior = info[var]['PP']

            varTR = info[var]['TR']
            varNF = info[var]['NF']
            varNR = info[var]['NR']

            FR.extend(freq)
            PP.extend(posterior)
            NR.extend(varNR)
            NF.extend(varNF)
            TR.extend(varTR)

        lineinfo['FR'] = FR
        lineinfo['PP'] = PP
        lineinfo['NF'] = NF
        lineinfo['NR'] = NR
        lineinfo['TR'] = TR

        qual = max([int(pp) for pp in lineinfo['PP']])
        linefilter = list(set(linefilter))
        vcfDataLine = {'chrom':chrom,'pos':POS,'ref':ref,'alt':alt,'id':id,'info':lineinfo,'filter':linefilter,'qual':qual,'format':format}

        # Filtering which can only be done at this level.
        if 'ma' in vcfFile.getfilter() and nVariants > 1:
            vcfDataLine['filter'].append('ma')

        # Now pull out the genotype informations for each sample, and compute genotype likelihoods and
        # Phred-based qualities scores.
        theVar = variants[0]
        testHap = call.haplotypes[0]
        refHap = Haplotype(testHap.refName, testHap.startPos, testHap.endPos, (), refFile, 0)
        varHap = Haplotype(testHap.refName, testHap.startPos, testHap.endPos, (theVar,), refFile, 0)

        AAGenotype = DiploidGenotype(2, [refHap,refHap])
        ABGenotype = DiploidGenotype(2, [refHap,varHap])
        BBGenotype = DiploidGenotype(2, [varHap,varHap])

        for i from 0 <= i < nIndividuals:

            genotype = genotypeCalls[i]

            # Missing genotype call. Probably due to zero coverage for this individual.
            if genotype is None:
                vcfDataLine[samplesThisPop[i]] = dict(GT=[[".", "/", "."]], GL=[-1,-1,-1], GQ=[-1], NR=[0])
                continue

            # Have genotype. Need to construct call string ("0/1" etc), and compute likelihoods
            # and quality.
            else:
                genotypeQuality = int((-10.0 * log10(max(1e-10, 1.0-genotypeLikelihoods[i]))))

                # Construct genotype string for un-phased genotypes
                hapOne = genotype.hap1
                hapTwo =  genotype.hap2

                varsHapOne = hapOne.variants
                varsHapTwo = hapTwo.variants

                indexOfCalledVarHapOne = 0
                indexOfCalledVarHapTwo = 0

                for varIndex from 0 <= varIndex < nVariants:

                    theVar = variants[varIndex]

                    if theVar in varsHapOne:
                        indexOfCalledVarHapOne = varIndex + 1 # 0 is the reference

                    if theVar in varsHapTwo:
                        indexOfCalledVarHapTwo = varIndex + 1 # 0 is the reference

                GT = [str(indexOfCalledVarHapOne), "/", str(indexOfCalledVarHapTwo)]

                # Genotype likelihoods.
                genotypeLikelihoodsThisInd = None

                if nVariants > 1:
                    genotypeLikelihoodsThisInd = [-1,-1,-1]
                else:
                    theBuffer = call.readBuffers[i]
                    nReadsThisIndividual = theBuffer.windowEnd - theBuffer.windowStart

                    AALike = round(AAGenotype.calculateDataLikelihood(theBuffer.windowStart, theBuffer.windowEnd, i, nReadsThisIndividual, nIndividuals), 2)
                    ABLike = round(ABGenotype.calculateDataLikelihood(theBuffer.windowStart, theBuffer.windowEnd, i, nReadsThisIndividual, nIndividuals), 2)
                    BBLike = round(BBGenotype.calculateDataLikelihood(theBuffer.windowStart, theBuffer.windowEnd, i, nReadsThisIndividual, nIndividuals), 2)

                    genotypeLikelihoodsThisInd = [AALike, ABLike, BBLike]

                    if options.verbosity >= 3:
                        logger.debug("Ref hap = %s" %(refHap))
                        logger.debug("Var hap = %s" %(varHap))
                        logger.debug("Genotype likelihoods for genotypes (\n%s, \n%s, \n%s) are as follows AA=%s, AB=%s, BB=%s" %(AAGenotype,ABGenotype,BBGenotype,AALike,ABLike,BBLike))

                vcfDataLine[samplesThisPop[i]] = dict(GT=[GT], GL=genotypeLikelihoodsThisInd, GQ=[genotypeQuality], NR=[nReadsThisIndividual])

        # Write variant call info to VCF file.
        vcfFile.write_data(outputFile, vcfDataLine)


cdef tuple refAndAlt(char* chrom, int POS, list variants, FastaFile refFile):
    """
    Calculates the right REF and ALT fields for a set of variants.
    In general, this is pretty difficult - it's ok if there's just
    one, but when there's more then one, it's a total nightmare.
    http://1000genomes.org/wiki/doku.php?do=show&id=1000_genomes%3Aanalysis%3Avcf4.0
    """
    cdef int nonSnpPresent = 0
    cdef int indelPresent = 0
    cdef Variant v
    cdef list ALT
    cdef list seq
    cdef bytes REF

    for v in variants:

        nAdded = v.nAdded
        nRemoved = v.nRemoved

        if nRemoved != 1 or nAdded != 1:
            nonSnpPresent = 1

            if nRemoved != nAdded:
                indelPresent = 1

    if not nonSnpPresent:

        # Just SNPs
        REF = refFile.getCharacter(chrom, POS)
        ALT = [ v.added for v in variants ]
        return REF,ALT

    else:

        # If there are indels around, then the reference base is the base before the insertion
        # or the base before the first deleted base. This means that the alternative allele for
        # snps will have two bases

        logger.debug("Outputting Variants at pos %s. Variants are %s" %(POS, variants))

        rlen = max([v.nRemoved for v in variants]) # number of extra bases in REF
        REF = None

        if indelPresent:
            REF = refFile.getSequence(chrom, POS, (POS+rlen) + 1)
        else:
            REF = refFile.getSequence(chrom, POS, (POS+rlen))

        ALT = []

        for v in variants:
            seq = list(REF)

            if v.nRemoved == v.nAdded:
                seq[0: len(v.added)] = v.added
            else:
                seq[1:1+v.nRemoved] = v.added

            ALT.append("".join(seq))

        return REF,ALT
