import os, sys, glob

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext as build_pyx

extModules = []
corMods = ['../../../coreutils/chaplotype.pxd', '../../../coreutils/variant.pxd', '../../../coreutils/fastafile.pxd', '../../../coreutils/calign.pxd', '../../../coreutils/samtoolsWrapper.pxd']
incDirs = ["../../../coreutils/samtools", "../../../coreutils", "./"]
incDirs.append("../cortex_var/include/basic")
incDirs.append("../cortex_var/include/basic/event_encoding/base_encoding")
#incDirs.append("../cortex_var/include/event_encoding/solid_colour_encoding")
incDirs.append("../cortex_var/include/hash_table")
incDirs.append("../cortex_var/include/hash_table/open_hash")
incDirs.append("../cortex_var/include/cortex_con")
incDirs.append("../cortex_var/include/cortex_var/many_colours")
incDirs.append("../cortex_var/include/cortex_var/core")
incDirs.append("../cortex_var/include/")

cFlags = ["-msse2", "-msse3", "-funroll-loops", "-D_LARGEFILE64_SOURCE", "-D_FILE_OFFSET_BITS=64"]
cortexCFlags = ["-DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=1", "-DNUMBER_OF_COLOURS=1"]

extModules.append(Extension(name='cgenotype', sources=corMods+['cgenotype.pyx',], include_dirs=incDirs,  extra_compile_args=cFlags))
extModules.append(Extension(name='cpopulation', sources=corMods+['cpopulation.pyx'],  include_dirs=incDirs, extra_compile_args=cFlags))
extModules.append(Extension(name='cwindow', sources=corMods+['cwindow.pyx'], include_dirs=incDirs,  extra_compile_args=cFlags))
extModules.append(Extension(name='cfilter', sources=corMods+['cfilter.pyx'], include_dirs=incDirs,  extra_compile_args=cFlags))
extModules.append(Extension(name='variantFilter', sources=corMods+['variantFilter.pyx'], include_dirs=incDirs,  extra_compile_args=cFlags))
extModules.append(Extension(name='variantcaller', sources=corMods+['variantcaller.pyx', 'cpopulation.pxd'], include_dirs=incDirs,  extra_compile_args=cFlags))

cortexLibDirs = ["../cortex_var/src/obj/cortex_var/many_colours/"]

cortexObjs = []
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/cmd_line.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/binary_kmer.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/cortex_var.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/db_differentiation.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/dB_graph.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/dB_graph_population.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/db_variants.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/element.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/event_encoding.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/file_reader.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/graph_info.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/hash_table.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/hash_value.o")
cortexObjs.append("../cortex_var/src/obj/cortex_var/many_colours/seq.o")

#samtoolsSources = [a for a in glob.glob(os.path.join("../../../coreutils/samtools/", "*.c")) if "bgzip" not in a]
#samtoolsObj = [a.split("/")[-1].replace(".c", ".o") for a in samtoolsSources]

#extModules.append(Extension(name='cortexWrapper', sources=list(set(corMods + samtoolsSources)) + ['cortexWrapper.pyx', 'cortex.c'], include_dirs=incDirs, library_dirs=cortexLibDirs, extra_objects=cortexObjs, extra_compile_args=cFlags + cortexCFlags, libraries=['z'], language='c'))

setup(name = "PlatypusCythonModules", ext_modules=extModules, cmdclass = {'build_ext': build_pyx})
