import cython
cimport cython

cimport chaplotype
cimport fastafile

from chaplotype cimport Haplotype
from fastafile cimport FastaFile

cdef list getFilteredHaplotypes(Haplotype refHaplotype, list variants, int nVars, FastaFile refFile, str windowChr, int windowStart, int windowEnd, int maxHaplotypes, int maxReadLength)
