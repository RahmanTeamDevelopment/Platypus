"""
Code for identifying variants in illumina reads, based on Gerton's haplotype realignment
algorithm and initial implementation.
"""

from __future__ import division

import logging
import bamfileutils
import multiprocessing
import variantcaller
import extendedoptparse
import os


logger = logging.getLogger("Log")


def chromAndPosSort(x, y):
    """
    Comparison function for use in sort routines. Compares strings of the form
    chr10:0-100. Sorting is done first by chromosome, in alphabetical order, and then
    by start position in numerical order.
    """
    xChrom = x.split(":")[0]
    yChrom = y.split(":")[0]
    xStart = int(x.split(":")[1].split("-")[0])
    yStart = int(y.split(":")[1].split("-")[0])

    return cmp(xChrom, yChrom) or cmp(xStart, yStart)


def printRegions(args):
    """
    Split the reference sequence into sensibly-sized regions, and
    print them out.
    """
    parser = extendedoptparse.OptionParser()

    parser.add_option("--refFile",dest="refFile", help="Fasta file of reference. Index must be in same directory", action='store', type='string', required=True)
    parser.add_option("--regions", dest="regions", type="list", help = "region as comma-separated list of chr:start-end, or just list of chr, or nothing", default=None, required=False, action = 'store')
    parser.add_option("--bamFiles", dest="bamFiles", type="list", help = "Comma-delimited list of bam file names", default=None, required=True)
    parser.add_option("--processRegionSize", dest="processRegionSize", type="int", help = "Max size of region to process. Chromosomes will be broken up into regions of this size", default=100000, required=False)
    parser.add_option("--parseNCBI", dest="parseNCBI", help="", type=int, action='store', default=0)

    (options, args) = parser.parse_args(args)
    regions = bamfileutils.getRegions(options)

    for region in regions:
        print "%s:%s-%s" %(region[0],region[1],region[2])


def runVariantCaller(options):
    """
    Run the indel caller

    We feed in a series of windows, each of which has a list of all the possible haplotypes in
    that window. For N P-ploid individuals, there are NP haplotypes to choose and we calculate
    the collection of NP haplotypes which maximise the likelihood.
    """
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    log = logging.getLogger('Log')

    ch = logging.StreamHandler()
    fh = logging.FileHandler(options.logFileName, 'w')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    if options.verbosity == 0:
        log.setLevel(logging.ERROR)
        ch.setLevel(logging.ERROR)
        fh.setLevel(logging.ERROR)
    elif options.verbosity == 1:
        log.setLevel(logging.WARNING)
        ch.setLevel(logging.WARNING)
        fh.setLevel(logging.WARNING)
    elif options.verbosity == 2:
        log.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
        fh.setLevel(logging.INFO)
    elif options.verbosity >= 3:
        # Debug goes to file only.
        log.setLevel(logging.DEBUG)
        ch.setLevel(logging.INFO)
        fh.setLevel(logging.DEBUG)
    else:
        raise StandardError, "Value of 'verbosity' input parameter must be between 0 and 3 inclusive"

    log.addHandler(ch)
    log.addHandler(fh)

    logger.info("Beginning variant calling")
    logger.info("Output will go to %s" %(options.output))

    fileNames = []

    if options.nCPU == 1:

        class DummyLock(object):
            def acquire(self):
                pass
            def release(self):
                pass

        # Don't multi-process, just use a single process. This is probably
        # more efficient, and certainly better for profiling.
        regions = bamfileutils.getRegions(options)

        # Write calls to temporary files, and then merge these into a
        # sorted vcf file at the end.
        lock = DummyLock()

        for region in regions:
            chromosome,start,end = region[::]
            fileName = options.output + "_temp_%s:%s-%s" %(chromosome, start, end)
            fileNames.append(fileName)
            variantcaller.generateGenotypesAndVariantsInAllRegions((region, fileName, options, lock))
    else:
        # Create process manager
        regions = bamfileutils.getRegions(options)
        manager = multiprocessing.Manager()
        lock = manager.Lock()
        theArgs = []

        for region in regions:
            chromosome,start,end = region[::]
            fileName = options.output + "_temp_%s:%s-%s" %(chromosome, start, end)
            fileNames.append(fileName)
            theArgs.append((region, fileName, options, lock))

        pool = multiprocessing.Pool(processes=options.nCPU)
        pool.map(variantcaller.generateGenotypesAndVariantsInAllRegions, theArgs)
        pool.close()
        pool.join()

    # Final output file
    outputVCF = open(options.output, 'wb')

    fileNamesForSorting = []

    for f in fileNames:
        regionThisFile = f.split("_temp_")[1]
        fileNamesForSorting.append((regionThisFile, f))

    for index,fileName in enumerate(sorted(fileNames, cmp=chromAndPosSort)):

        inputVCF = open(fileName, 'rb')

        for line in inputVCF:

            if line[0] == "#" and index == 0:
                outputVCF.write(line)
            elif line[0] == "#":
                continue
            else:
                outputVCF.write(line)

        # Close and remove temporary files
        inputVCF.close()

    for fileName in fileNames:
        os.remove(fileName)

    outputVCF.close()
    logger.info("Finished variant calling")

###################################################################################################

def callVariants(args):
    """
    Run the Platypus variant-caller, with the specified arguments
    """
    parser = extendedoptparse.OptionParser()

    parser.add_option("-o", "--output", dest="output",  help="Output SNP data file", action='store', type='string', default="AllVariants.vcf")
    parser.add_option("-n", "--nIndividuals", dest="nInd", help="Number of Individuals sequenced", action='store', type='int', default=1)
    parser.add_option("-p", "--ploidy", dest="ploidy", help="Number of copief of each chromosome per individual", action='store', type='int', default=2)
    parser.add_option("--refFile",dest="refFile", help="Fasta file of reference. Index must be in same directory", action='store', type='string', required=True)
    parser.add_option("--regions", dest="regions", type="list", help = "region as comma-separated list of chr:start-end, or just list of chr, or nothing", default=None, action = 'store')
    parser.add_option("--bamFiles", dest="bamFiles", type="list", help = "Comma-delimited list of bam file names", default=None, required=True)
    parser.add_option("--processRegionSize", dest="processRegionSize", type="int", help = "Size of region to use in each process. Chromosomes will be broken up into regions of this size and run in separate processes (or consecutively, if nCPU == 1)", default=30000000, required=False)
    parser.add_option("--bufferSize", dest="bufferSize", type="int", help = "Data will be buffered in regions of this size", default=1000000, required=False)
    parser.add_option("--minReads", dest="minReads", help="Minimum required number of reads in window", action='store', type='int', default=2)
    parser.add_option("--maxHaplotypes", dest="maxHaplotypes", help="Maximium haplotypes to consider in a given window", action='store', type='int', default=256)
    parser.add_option("--maxVariants", dest="maxVariants", help="Maximium variants to consider in a given window", action='store', type='int', default=8)
    parser.add_option("--maxReads", dest="maxReads", help="Maximium coverage in window", action='store', type='float', default=5000000)
    parser.add_option("--maxSize", dest="maxSize", help="Largest indel to consider", action='store', type='int', default=250)
    parser.add_option("--verbosity", dest="verbosity", help="Level of logging", action='store', type='int', default=2)
    parser.add_option("--minMapQual", dest="minMapQual", help="Minimum mapping quality of read. Any reads with map qual below this are ignored", action='store', type = 'int', default=20, required=False)
    parser.add_option("--maxBadQualBases", dest="maxBadQualBases", help="Max number of bases per read that can have base-quality <= 20", action='store', type = 'int', default=20, required=False)
    parser.add_option("--minGoodQualBases", dest="minGoodQualBases", help="Minimum number of bases per read that must have base-quality >= 20", action='store', type = 'int', default=20, required=False)
    parser.add_option("--minBaseQual", dest="minBaseQual", help="Minimum allowed base-calling quality. Any bases with qual below this are ignored in SNP-calling", action='store', type = 'int', default=20, required=False)
    parser.add_option("--maxReadLength", dest="rlen", help="Maximum read length", action='store', type = 'int', default=150)
    parser.add_option("--minFlank", dest="minFlank", help="Flank size for indel candidates", action='store', type = 'int', default=3)
    parser.add_option("--logFileName", dest="logFileName", help="Name of log file", action='store', type='string', default="log.txt")
    parser.add_option("--labels", dest="labels", help="regex to extract names from filenames", action='store', type = 'string', default= None)
    parser.add_option("--source", dest="sourceFile", help="vcf file(s) to get candidates from", action='store', type='list', default= None)
    parser.add_option("--freqAsPrior", dest="freqAsPrior", help="If 1, use the frequency of input variants as a prior", action='store', type = 'int', default=0)
    parser.add_option("--dataType", dest="dataType", help="Are we looking at individual,trio,population, or pool data?", required=False, action='store', default="population", type='choice', choices=["population"])
    parser.add_option("--nCPU", dest="nCPU", help="Number of processors to use", action='store', type='int', default=1)
    parser.add_option("--getVariantsFromBAMs", dest="getVariantsFromBAMs", help="If set to TRUE (default), variant candidates will be generated from BAMs as well as any other inputs", action='store', type='int', default=1)
    parser.add_option("--genSNPs", dest="genSNPs", help="If set to TRUE (default), SNP candidates will be considered", action='store', type='int', default=1)
    parser.add_option("--genIndels", dest="genIndels", help="If set to TRUE (default), Indel candidates will be considered", action='store', type='int', default=1)
    parser.add_option("--parseNCBI", dest="parseNCBI", help="", type=int, action='store', default=0)
    parser.add_option("--callOnlyIndels", dest="callOnlyIndels", help="If set to TRUE (default), only windows containing Indel candidates will be considered", action='store', type='int', default=0)
    parser.add_option("--strandFilter", dest="strandFilter", help="If set to TRUE only variants occuring at least once on each strand will be considered.", action='store', type='int', default=0)
    parser.add_option("--minPosterior", dest="minPosterior", help="Only variants with posterior >= this will be outpu to the VCF. Value is a Phred-score.", action='store', type='int', default=5)
    parser.add_option("--sbThreshold", dest="sbThreshold", help="P-value for strand-bias filtering..", action='store', type='float', default=1e-2)
    parser.add_option("--abThreshold", dest="abThreshold", help="P-value for allele-bias filtering..", action='store', type='float', default=1e-3)
    parser.add_option("--printVarsAndExit", dest="printVarsAndExit", help="If 1, print a list of variant candidates, and exit.", action='store', type='int', default=0)
    parser.add_option("--useCortex", dest="useCortex", help="If 1, Cortex will be used to assemble variant candidates for Platypus to call.", action='store', type='int', default=0)
    parser.add_option("--badReadsWindow", dest="badReadsWindow", help="Size of window around variant to look for low-quality bases.", action='store', type='int', default=11)
    parser.add_option("--badReadsThreshold", dest="badReadsThreshold", help="Variants where the median minimum quality in a window of badReadsWindow around the variant position falls below this value will be filtered with the flag 'badReads'.", action='store', type='int', default=15)

    (options, args) = parser.parse_args(args)
    runVariantCaller(options)

###################################################################################################
