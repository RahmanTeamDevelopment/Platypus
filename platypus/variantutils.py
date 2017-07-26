from __future__ import division

import logging
import pysam
from platypus.variant import Variant


logger = logging.getLogger("Log")


class VariantCandidateReader(object):
    """
    A class to read variant information from a vcf file, and return
    a stream of variant instances.
    """
    def __init__(self, fileNames):
        """
        Constructor. Takes the name of a vcf file.
        """
        for fileName in fileNames:
            if ".gz" not in fileName:
                try:
                    pysam.tabix_compress(fileName, fileName + ".gz")
                    fileName += ".gz"
                except IOError:
                    pass

            try:
                pysam.tabix_index(fileName, preset="vcf")
            except IOError:
                pass

            self.vcfFiles = [pysam.Tabixfile(fileName) for fileName in fileNames]

    def __del__(self):
        """
        Destructor. Make sure to close files.
        """
        pass

    def Variants(self, chromosome, start, end, freqAsPrior=False):
        """
        Generator funtion. Yields variants in order of
        genomic co-ordinate.
        """
        varList = []

        for vcfFile in self.vcfFiles:
            try:
                vcfLines = vcfFile.fetch(chromosome, start, end, parser=pysam.asVCF())
            except Exception, e:
                logger.warning("Could not retrieve variants from source file in region %s:%s-%s. Error was %s" %(chromosome,start,end,e.message))
                continue

            for line in vcfLines:

                pos = line.pos
                ref = line.ref
                alt = line.alt
                info = line.info

                infoEntries = info.split(";")
                freq = None

                if freqAsPrior:
                    for entry in infoEntries:

                        cols = entry.split("=")

                        if cols[0] == "AF":
                            if len(cols) == 2:
                                freq = float(cols[1])
                            else:
                                # Skip multi-allelic sites for now
                                continue

                lenRef = len(ref)
                lenAlt = len(alt)

                # SNP
                if lenRef == 1 and lenAlt == 1:
                    var = Variant(chromosome, pos, ref, alt, 0, 0, 0, 0, 0)

                    if freqAsPrior:
                        var.freqPrior = freq

                    varList.append(var)

                # Anything else
                else:

                    # VCF4 is -1 indexed for indels, so trim off first base
                    ref = ref[1:]
                    alt = alt[1:]
                    removed = ref
                    added = alt

                    # Trim the matching bits off and shift position. This will decompose
                    # multi-variant sites into individual alleles at different positions.
                    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
                        ref = ref[1:]
                        alt = alt[1:]
                        removed = ref
                        added = alt
                        pos +=1

                    var = Variant(chromosome, pos, removed, added, 0, 0, 0, 0, 0)

                    if freqAsPrior:
                        var.freqPrior = freq

                    varList.append(var)

        varList = sorted(list(set(varList)))
        logger.info("Found %s variants in region %s in source file" %(len(varList), "%s:%s-%s" %(chromosome,start,end)))
        return varList

###################################################################################################
