import glob
import os

from distutils.command.build_clib import build_clib
from setuptools import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize


ext_modules=[
    Extension("demo", ["demo.pyx"])
]

include_dirs = [
    "platypus",
    "platypus/c",
    "platypus/c/samtools"
]

cython_directives = {
    "boundscheck": False,
    "nonecheck" : False,
    "cdivision" : True,
    "profile" : False,
    "initializedcheck" : False,
    "wraparound" : True
}

samtools_source_files = glob.glob(
    os.path.join("platypus", "c", "samtools", "*.c")
)

optimisation_args = [
    "-msse2",
    "-msse3",
    "-funroll-loops"
]

large_file_opts = [
    "-D_LARGEFILE64_SOURCE",
    "-D_FILE_OFFSET_BITS=64"
]

compile_args = optimisation_args + large_file_opts

libsamtools = (
    'samtools', {
        'sources': samtools_source_files
    }
)


modules = [
    Extension(
        name='platypus.samtoolsWrapper',
        sources=['platypus/samtoolsWrapper.pyx', 'platypus/c/pysam_util.c'],
        include_dirs=include_dirs,
        libraries=['z'],
        language='c',
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.fastafile',
        sources=['platypus/fastafile.pyx'],
        include_dirs=include_dirs,
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.bamfileutils',
        sources=['platypus/bamfileutils.pyx'],
        include_dirs=include_dirs,
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.variant',
        sources=['platypus/variant.pyx'],
        include_dirs=include_dirs,
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.cerrormodel',
        sources=['platypus/cerrormodel.pyx', 'platypus/c/tandem.c'],
        include_dirs=include_dirs,
        language='c',
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.calign',
        sources=['platypus/calign.pyx', 'platypus/c/align.c'],
        include_dirs=include_dirs,
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.chaplotype',
        sources=['platypus/chaplotype.pyx'],
        include_dirs=include_dirs,
        language='c',
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.cgenotype',
        sources=['platypus/cgenotype.pyx'],
        include_dirs=include_dirs,
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.cpopulation',
        sources=['platypus/cpopulation.pyx'],
        include_dirs=include_dirs,
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.cwindow',
        sources=['platypus/cwindow.pyx'],
        include_dirs=include_dirs,
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.cfilter',
        sources=['platypus/cfilter.pyx'],
        include_dirs=include_dirs,
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.variantFilter',
        sources=['platypus/variantFilter.pyx'],
        include_dirs=include_dirs,
        extra_compile_args=compile_args
    ),
    Extension(
        name='platypus.variantcaller',
        sources=['platypus/variantcaller.pyx'],
        include_dirs=include_dirs,
        extra_compile_args=compile_args
    )
]

setup(
    name="Platypus",
    version='0.1.5',
    description="Small variant calling",
    url='https://github.com/RahmanTeamDevelopment/Platypus',
    author='RahmanTeam',
    author_email='rahmanlab@icr.ac.uk',
    license='MIT',
    libraries=[libsamtools],
    cmdclass={'build_clib': build_clib, 'build_ext': build_ext},
    ext_modules = cythonize(
        modules,
        compiler_directives=cython_directives
    ),
    packages=[
        'platypus',
    ],
    scripts=[
        "bin/Platypus.py",
        "bin/platypus",
        "test/smoke/check_installation_succeeded.bash",
    ],
    zip_safe=False,
    include_package_data=True,
    install_requires=[
        "Cython==0.25.2",
        "pysam==0.10.0"
    ]
)