
from platypus.chaplotype cimport Haplotype
from platypus.fastafile cimport FastaFile

cdef list getFilteredHaplotypes(Haplotype refHaplotype, list variants, int nVars, FastaFile refFile, str windowChr, int windowStart, int windowEnd, int maxHaplotypes, int maxReadLength)
