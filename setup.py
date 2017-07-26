import os, sys, glob

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext as build_pyx

extModules = []

extModules.append(Extension(name='samtoolsWrapper', sources=['samtoolsWrapper.pyx', 'pysam_util.c'] + glob.glob(os.path.join("samtools", "*.c")), include_dirs=["samtools", "./"],libraries=['z'], language='c', extra_compile_args=["-D_LARGEFILE64_SOURCE", "-D_FILE_OFFSET_BITS=64", "-msse2", "-msse3", "-funroll-loops"]))

extModules.append(Extension(name='fastafile', sources=['fastafile.pyx'], include_dirs=["samtools", "./"],  extra_compile_args=["-msse2", "-msse3", "-funroll-loops"]))

extModules.append(Extension(name='bamfileutils', sources=['bamfileutils.pyx', 'samtoolsWrapper.pxd'], include_dirs=["samtools", "./"], extra_compile_args=["-msse2", "-msse3", "-funroll-loops"]))

extModules.append(Extension(name='variant', sources=['variant.pyx', 'samtoolsWrapper.pxd'], include_dirs=["samtools", "./"], extra_compile_args=["-msse2", "-msse3", "-funroll-loops"]))

extModules.append(Extension(name='histogram', sources=['histogram.pyx'], include_dirs=["samtools", "./"], extra_compile_args=["-msse2", "-msse3", "-funroll-loops"]))

extModules.append(Extension(name='cerrormodel', sources=['cerrormodel.pyx', 'tandem.c'], include_dirs=["./"], language='c', extra_compile_args=["-msse2", "-msse3", "-funroll-loops"]))

extModules.append(Extension(name='calign', sources=['calign.pyx', 'align.c', 'samtoolsWrapper.pxd', 'cerrormodel.pxd'], include_dirs=["samtools", "./"], extra_compile_args=["-msse2", "-msse3", "-funroll-loops"]))

extModules.append(Extension(name='chaplotype', sources=['chaplotype.pyx', 'variant.pxd','samtoolsWrapper.pxd', 'calign.pxd'], include_dirs=["samtools", "./"], language='c', extra_compile_args=["-msse2", "-msse3", "-funroll-loops"]))


# Set-up and install python modules
setup(name="CoreUtilsCythonModules", ext_modules=extModules, cmdclass={'build_ext': build_pyx})
setup(name="filez", py_modules=['filez'])
setup(name="extendedoptparse", py_modules=['extendedoptparse'])



from distutils.core import setup

setup(name="window", py_modules=['window'])
setup(name="variantutils", py_modules=['variantutils'])
setup(name="platypusexceptions", py_modules=['platypusexceptions'])
setup(name="extendedoptparse", py_modules=['extendedoptparse'])
setup(name="runner", py_modules=['runner'])