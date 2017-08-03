
************
Introduction
************

Platypus is a tool for identifying small variants (SNPs and indels) in BAM files. This code is a fork of Platypus version 0.1.5, which is used
in OpEx. The latest Platypus code is available on `Github <https://github.com/andyrimmer/Platypus>`_.


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

.. csv-table::
	:header: "Option", "Type", "Default Value", "Example Usage", "Description"
	:delim: |

	-h --help |  NA | NA | `--help` | If this option is set then Platypus will print a help message and then exit
	--bamFiles |  Comma-separated list of strings | None | `--bamFiles=test1.bam,test2.bam` | A list of BAM file names. All the specified BAM files will be searched for variants. The output will go into a single VCF file
	--refFile | String | None | `--refFile=human_37.fa` | Name of the input reference FASTA file. This file must be indexed (using `Samtools faidx` for example) and must contain the same reference sequence that was used to map the reads in the input BAM file(s)


More text goes here.

