
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

	--processRegionSize | Integer, 30,000,000 | `--processRegionSize=10000000` | Platypus breaks up the genome into regions of this size and processes them in parallel (if nCPU > 1) or consecutively (if nCPU == 1)
	--bufferSize | Integer | 1,000,000 | `--bufferSize=100000` | The maximum size of region (per process, and as a genomic interval) that Platypus will read into memory at any point. This can be used to control memory usage
	--maxReads | Integer | 5,000,000 | `--maxReads=100000` | This sets an upper limit on the amount of data that Platypus will try to process in one region. If the number of reads in a region of `bufferSize` is larger than this then Platypus will skip the region. This option can be used, in conjunction with `bufferSize` to control memory usage.
	--maxBadQualBases | Integer | ?  | ?

--maxBadQualBases=MAXBADQUALBASES
                      Max number of bases per read that can have base-
                      quality <= 20
--minGoodQualBases=MINGOODQUALBASES
                      Minimum number of bases per read that must have base-
                      quality >= 20
--maxReadLength=RLEN  Maximum read length
--labels=LABELS       regex to extract names from filenames

Variant calling options
=======================

.. csv-table::
	:header: "Option", "Type", "Default Value", "Example Usage", "Description"
	:delim: |

	--genSNPs | Boolean | 1 (True) | `--genSNPs=1` | If set to 1, Platypus will call SNPs. If set to 0 Platypus will not call SNPs
	--genIndels | Boolean | 1 (True) | `--genIndels=1` | If set to 1, Platypus will call Indels. If set to 0 Platypus will not call Indels
	--minBaseQual | Integer | 20 | `--minBaseQual=25` | Only bases with base quality >= this value will be examined when generating the initial list of SNP candidates
	--minMapQual | Integer | 20 | `--minMapQual=25` | Only bases with base quality >= this value will be examined when generating the initial list of SNP candidates

--maxHaplotypes=MAXHAPLOTYPES
                      Maximium haplotypes to consider in a given window
--maxVariants=MAXVARIANTS
                      Maximium variants to consider in a given window
--maxSize=MAXSIZE     Largest indel to consider
--minReads=MINREADS   Minimum required number of reads in window
--callOnlyIndels=CALLONLYINDELS
                      If set to TRUE (default), only windows containing
                      Indel candidates will be considered
--strandFilter=STRANDFILTER
                      If set to TRUE only variants occuring at least once on
                      each strand will be considered.
--getVariantsFromBAMs=GETVARIANTSFROMBAMS
                      If set to TRUE (default), variant candidates will be
                      generated from BAMs as well as any other inputs

--minFlank=MINFLANK   Flank size for indel candidates


Variant filtering options
=========================

.. csv-table::
	:header: "Option", "Type", "Default Value", "Example Usage", "Description"
	:delim: |

--minPosterior=MINPOSTERIOR
                      Only variants with posterior >= this will be outpu to
                      the VCF. Value is a Phred-score.
--sbThreshold=SBTHRESHOLD
                      P-value for strand-bias filtering..
--abThreshold=ABTHRESHOLD
                      P-value for allele-bias filtering..
--badReadsWindow=BADREADSWINDOW
                      Size of window around variant to look for low-quality
                      bases.
--badReadsThreshold=BADREADSTHRESHOLD
                      Variants where the median minimum quality in a window
                      of badReadsWindow around the variant position falls
                      below this value will be filtered with the flag
                      'badReads'


Miscellaneous options
=====================

.. csv-table::
	:header: "Option", "Type", "Default Value", "Example Usage", "Description"
	:delim: |

	--source | String | None | `--source=thousand_genomes_snps.vcf.gz` | Name of an input VCF file to be used as a source of variant candidates. See :ref:`Supplying variant candidates from VCF<supplying_variant_candidates_from_vcf>` for more details

--source=SOURCEFILE   vcf file(s) to get candidates from
--freqAsPrior=FREQASPRIOR
                      If 1, use the frequency of input variants as a prior
--parseNCBI=PARSENCBI
--printVarsAndExit=PRINTVARSANDEXIT
                      If 1, print a list of variant candidates, and exit.


Deprecated command-line options
===============================

The following command-line options are deprecated, and should not be used.

.. csv-table::
	:header: "Option", "Description"
	:delim: |

	-n --nIndividuals | Was used to set the number of individuals in the input BAM files. Now each BAM file is assumed to contain data from only one individual
	-p --ploidy | Was used to set the ploidy of the samples in the BAM files. Now this is fixed at 2.
	--dataType | Was used to distinguish between individual, trio and pooled sequencing datasets


.. _specifying_calling_regions:

Specifying calling regions
==========================

It often useful to call variants on only a subset of the genome, e.g. a particular gene or exon. Platypus supports this mode through the `--regions` option. This option can be used in several ways:

	* A comma-separated list of chromosome:start-end coordinates can be given, e.g. `--regions=chr1:0-100,chr2:300-400,chrX:1000-20000`. Platypus will search for variants in just these regions
	* A comma-separated list of chromosomes can be given, e.g. `--regions=chr1,chr2,chr20,chrY`. Platypus will search for variants in just these chromosomes, and the lengths of the chromosomes will be taken from either the BAM file header or the reference FASTA file index
	* A text file name may be specified, e.g. `--regions=my_regions.txt`. This file must contain only lines with the format `chrom:start-end`, e.g. `chr1:1000-20000` with one region per line.
	* If no regions are specified then Platypus will call variants across the whole genome. It will check the BAM header file and the reference FASTA file to determine the list of chromosomes and other contigs, as well as their lengths


.. _supplying_variant_candidates_from_vcf:

Supplying variant candidates from VCF
=====================================
