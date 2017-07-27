
************
Introduction
************

Platypus is a tool for identifying small variants (SNPs and indels) in BAM files.


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

To install and run Platypus v1.3.0 you will need `Git <https://git-scm.com>`_, `Python <https://www.python.org>`_ (only
version 2.7.X is currently supported), and `Virtualenv <https://virtualenv.pypa.io/en/stable/>`_. 


*****************
Running Platypus
*****************

Once downloaded and correctly installed, Platypus can be run with the following simple command::

    env/bin/platypus --bamFiles=input.bam --refFile=reference.fa --output=variant_calls.vcf

By default, Platypus takes three command line arguments: the name of the input BAM file, the of input reference file and the name of the output VCF file. 

* The input bam file (-i) should follow the `BAM format <http://samtools.github.io/hts-specs/SAMv1.pdf>`_ and should contain the mapped reads. The .bai index file should also be present in the same directory.
