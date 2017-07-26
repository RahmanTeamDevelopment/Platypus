cimport fastafile
from fastafile cimport FastaFile

cdef list filterVariants(list varList, FastaFile refFile, int maxReadLength, int minSupport, int strandFilter, int maxDiff, int verbosity)
