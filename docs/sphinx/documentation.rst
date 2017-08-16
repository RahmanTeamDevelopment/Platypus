
************
Introduction
************

Platypus is a tool for identifying small variants (SNPs and indels) in BAM files. This code is a fork of Platypus version 0.1.5, which is used
in OpEx. If you are looking for the most recent version of Platypus, the code is available `here <https://github.com/andyrimmer/Platypus>`_.


************
Installation
************

Stable releases of Platypus can be downloaded from Github `here <https://github.com/RahmanTeamDevelopment/Platypus/releases>`_
in either ``.zip`` or ``.tar.gz`` format. To unpack these run one of the following commands::

    unzip Platypus.X.X.zip

or::

    tar -xvzf Platypus.X.X.tar.gz

and then you can install Platypus with the following commands::

    cd Platypus
    ./install.sh

Platypus uses ``virtualenv`` and ``pip`` to manage all its extra dependencies, which means that it will not clutter up your system by installing
things globally. Everthing it installs will go into a sub-directory in the ``Platypus`` directory (specifically, ``Platypus_X.X.X/env``). If
you delete Platypus then everything it has installed will also be deleted. Once the installation script has finished successfully,
Platypus is ready for use. 


Dependencies
============

To install and run Platypus v0.1.5 you will need `Python <https://www.python.org>`_ (only version 2.7.X is currently supported),
and `Virtualenv <https://virtualenv.pypa.io/en/stable/>`_. 


*****************
Running Platypus
*****************

Once downloaded and correctly installed, Platypus can be run with the following simple command::

    env/bin/platypus callVariants --bamFiles=input.bam --refFile=reference.fa --output=variant_calls.vcf

This is the simplest way to run Platypus, specifying only the names of the input BAM file and input reference file. Variant
calls will then be outupt to a file called `AllVariants.vcf`. 

* The input bam file should follow the `BAM format <http://samtools.github.io/hts-specs/SAMv1.pdf>`_ and should contain mapped reads sorted by genomic coordinate. Platypus has been tested using BAM files produced by `Stampy <http://www.well.ox.ac.uk/project-stampy>`_ and `BWA <http://bio-bwa.sourceforge.net>`_, but should also work with other short read aligners. The BAM index file (.bai) should also be present in the same directory; this can be created using `Samtools <http://www.htslib.org>`_.


********************
Command Line Options
********************

Platypus has a large number of configuration options that can be specified on the command-line. Most of these can be ignored, as they take default values
and are already configured for optimal results. Below is a full descripton of all configuration options.


Basic options
=============

.. csv-table::
    :header: "Option", "Type", "Default Value", "Example Usage", "Description"
    :delim: |

    -h --help |  NA | NA | `--help` | If this option is set then Platypus will print a help message and then exit
    --bamFiles |  Comma-separated list of strings | None | `--bamFiles=test1.bam,test2.bam` | A list of BAM file names. All the specified BAM files will be searched for variants. The output will go into a single VCF file
    --refFile | String | None | `--refFile=human_37.fa` | Name of the input reference FASTA file. This file must be indexed (using `Samtools faidx` for example) and must contain the same reference sequence that was used to map the reads in the input BAM file(s)
    -o --output | String | AllVariants.vcf | `--output=variant_calls.vcf` | Name of the output VCF file
    --regions | Comma-separated list of strings or the name of a text file | None | `--regions=chr1:1000-2000,chr2:2000-4000` | Platypus will produce variant calls only in the specified regions. See :ref:`Specifying calling regions<specifying_calling_regions>` for a more detailed description 
    --verbosity | Integer from 0 to 3 | 2 | `--verbosity=3` | Sets the level of output logging for Platypus. Increase this value for more verbose log output
    --logFileName | String | log.txt | `--logFileName=platypus_log_file.txt` | Name of the output log file.
    --nCPU | Integer | 1 | `--nCPU=3` | Number of processes to use for variant calling. If > 1 then Platypus will run in multiple processes and merge the VCF files into one file at the end


.. _bam_data_quality_filtering:

BAM Data processing and filtering options
=========================================

Platypus does not necessarily use every read, in the BAM for variant calling. Whole reads may be filtered out based on various quality criteria, and segments of individual reads may also be filtered. The table below gives a summary of the configuration options that may be used to control read-level filtering.

.. csv-table::
    :header: "Option", "Type", "Default Value", "Example Usage", "Description"
    :delim: |

    --processRegionSize | Integer | 30,000,000 | `--processRegionSize=10000000` | Platypus breaks up the genome into regions of this size and processes them in parallel (if nCPU > 1) or consecutively (if nCPU == 1)
    --bufferSize | Integer | 1,000,000 | `--bufferSize=100000` | The maximum size of region (per process, and as a genomic interval) that Platypus will read into memory at any point. This can be used to control memory usage
    --maxReads | Integer | 5,000,000 | `--maxReads=100000` | This sets an upper limit on the amount of data that Platypus will try to process in one region. If the number of reads in a region of `bufferSize` is larger than this then Platypus will skip the region. This option can be used, in conjunction with `bufferSize` to control memory usage.
    --maxBadQualBases | Integer | 20  |  `--maxBadQualBases=10` | If a read has > this number of bases with quality < 20 then it will be filtered out and not used for variant calling.
    --minGoodQualBases | Integer | 20  |  `--minGoodQualBases=10` | If a read has < this number of bases with quality >= 20 then it will be filtered out and not used for variant calling.



Variant calling options
=======================

.. csv-table::
    :header: "Option", "Type", "Default Value", "Example Usage", "Description"
    :delim: |

    --genSNPs | Boolean | 1 (True) | `--genSNPs=1` | If set to 1, Platypus will call SNPs. If set to 0 Platypus will not call SNPs
    --genIndels | Boolean | 1 (True) | `--genIndels=1` | If set to 1, Platypus will call Indels. If set to 0 Platypus will not call Indels
    --minBaseQual | Integer | 20 | `--minBaseQual=25` | Only bases with base quality >= this value will be examined when generating the initial list of SNP candidates
    --minMapQual | Integer | 20 | `--minMapQual=25` | Only bases with base quality >= this value will be examined when generating the initial list of SNP candidates
    --maxVariants | Integer | 8 | `--maxVariants=12` | A limit on the number of variant candidates that Platypus will consider in any single window. For more detals see :ref:`Filtering variant candidates<filtering_variant_candidates>` for more details
    --maxHaplotypes | Integer | 256 | `--maxHaplotypes=500` | A limit on the number of haplotpyes that Platypus will consider in any single window. For more details see :ref:`Generating haplotypes<generating_haplotypes>` for more details
    --minReads | Integer | 1 | `--minReads=10` | There must be >= minReads reads with mapping quality >= minMapQual supporting any variant for it to be considered. For SNPs the base qualities must be >= minBaseQual 
    --maxSize | Integer | 250 | `--maxSize=250` | An upper limit on the size of indels that Platypus will call
    --getVariantsFromBAMs | Boolean | 1 | `--getVariantsFromBAMs=1` | If set to 0 then Playpus will not scan the BAM for variant candidates. This only makes sense if there are candidates coming from elsewhere, e.g. an input VCF of variants
    --strandFilter | Boolean | 0 | `--strandFilter=1` | If set to 1 then Platypus will require at least 1 read on each of the forward and reverse strands in order to consider a variant candidate
    --minFlank | Integer | 3 | `--minFlank=10` | Any indel candidate in the BAM must have >= minFlank reference-matching bases on either side of it to be considered



Variant filtering options
=========================

The VCF output of Platypus contains a number of metrics and flags that can be used to filter the variant calls. The following options can be used to control the filters.

.. csv-table::
    :header: "Option", "Type", "Default Value", "Example Usage", "Description"
    :delim: |

    --minPosterior | Float | 5 | `--minPosterior=10` | Only variant calls wth a posterior / QUAL value >= this will be output to the VCF. This value is a :ref:`Phred score<phred_scores>`.
    --sbThreshold | Float | 0.01 | `--sbThreshold=0.02` | P-value threshold for the strand-bias filter. Variants with a strand-bias p-value < this threshold will be flagged as `sb` in the VCF filter column. See :ref:`Filtering variant candidatess<filtering_variant_candidates>` for more details.
    --abThreshold | Float | 0.001 | `--abThreshold=0.002` | P-value threshold for the allele-bias filter. Variants with an allele-bias p-value < this threshold will be flagged as `ab` in the VCF filter column. See :ref:`Filtering variant candidates<filtering_variant_candidates>` for more details.
    --minPosterior | Float | 5 | `--minPosterior=10` | Only variant calls wth a posterior / QUAL value >= this will be output to the VCF. This value is a :ref:`Phred score<phred_scores>`.
    --badReadsWindow | Integer | 11 | `--badReadsWindow=7` | Size of window around variant location to check for low quality bases. If the median minimum quality in this window is < badReadsThreshold then the variant will be filtered with the flag 'badReads' in the output VCF. See :ref:`Filtering variant candidates<filtering_variant_candidates>` for more details.
    --badReadsThreshold | Integer | 15 | `--badReadsThreshold=20` | Threshold for the bad reads filter. See :ref:`Filtering variant candidates<filtering_variant_candidates>` for more details.



Miscellaneous options
=====================

.. csv-table::
    :header: "Option", "Type", "Default Value", "Example Usage", "Description"
    :delim: |

    --source | String | None | `--source=thousand_genomes_snps.vcf.gz` | Name of an input VCF file to be used as a source of variant candidates. See :ref:`Supplying variant candidates from VCF<supplying_variant_candidates_from_vcf>` for more details
    --freqAsPrior | Boolean | 0 | `--freqAsPrior=1` | Must be used in conjunction with `--source`. If this is set to 1 then Platypus will use the AF tag in the input VCF as a prior on the variant. See :ref:Supplying variant candidates from VCF<supplying_variant_candidates_from_vcf>` for more details
    --parseNCBI | Boolean | 0 | `--parseNCBI=1` | This must be set to 1 if the reference file has sequence names in the NCBI format
    --printVarsAndExit | Boolean | 0 | `--printVarsAndExit=1` | If set to 1, Platypus will simply print a list of all variant candidates in the calling regions, and exit without making calls


Deprecated command-line options
===============================

The following command-line options are deprecated, and should not be used.

.. csv-table::
    :header: "Option", "Description"
    :delim: |

    -n --nIndividuals | Was used to set the number of individuals in the input BAM files. Now each BAM file is assumed to contain data from only one individual
    -p --ploidy | Was used to set the ploidy of the samples in the BAM files. Now this is fixed at 2.
    --dataType | Was used to distinguish between individual, trio and pooled sequencing datasets
    --maxReadLength | Was used to limit the size of reads that Platypus would process
    --labels | Was used to specify how to convert BAM file names into the VCF sample column headers
    --callOnlyIndels | If set to 1 then Platypus will only call variants in regions where there are indel candidates


.. _specifying_calling_regions:

Specifying calling regions
==========================

It often useful to call variants on only a subset of the genome, e.g. a particular gene or exon. Platypus supports this mode through the `--regions` option. This option can be used in several ways:

    * A comma-separated list of chromosome:start-end coordinates can be given, e.g. `--regions=chr1:0-100,chr2:300-400,chrX:1000-20000`. Platypus will search for variants in just these regions
    * A comma-separated list of chromosomes can be given, e.g. `--regions=chr1,chr2,chr20,chrY`. Platypus will search for variants in just these chromosomes, and the lengths of the chromosomes will be taken from either the BAM file header or the reference FASTA file index
    * A text file name may be specified, e.g. `--regions=my_regions.txt`. This file must contain only lines with the format `chrom:start-end`, e.g. `chr1:1000-20000` with one region per line.
    * If no regions are specified then Platypus will call variants across the whole genome. It will check the BAM header file and the reference FASTA file to determine the list of chromosomes and other contigs, as well as their lengths


========================================
Description of Variant Calling Algorithm
========================================

Blah


.. _generating_variants:

Generating variants
===================

Description of how Platypus generates variant candidates from BAMs, input VCFs etc



.. _filtering_variant_candidates:

Filtering variant candidates
============================

Blah



.. _generating_haplotypes:

Generating haplotypes
=====================

Blah



.. _supplying_variant_candidates_from_vcf:

Supplying variant candidates from VCF
=====================================

Blah



.. _phred_scores:

Phred scores
============

Blah
