#cython: boundscheck=False
#cython: cdivision=True
#cython: nonecheck=False

"""
Cython module used to interface to Cortex routines.
"""

from __future__ import division

cimport cython

import logging
import samtoolsWrapper
cimport samtoolsWrapper
cimport variant

from samtoolsWrapper cimport cAlignedRead
from variant cimport Variant

logger = logging.getLogger("Log")

###################################################################################################

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *memset(void *buffer, int ch, size_t count)

###################################################################################################

cdef extern from "limits.h":
    cdef int LINE_MAX

###################################################################################################

cdef extern from "math.h":
    double exp(double)
    double log(double)

###################################################################################################

cdef extern from "binary_kmer.h":

    ctypedef struct BinaryKmer:
        pass

    ctypedef BinaryKmer* Key

    ctypedef struct KmerSlidingWindow:
        int nkmers
        BinaryKmer* kmer

    ctypedef struct KmerSlidingWindowSet:
        int nwindows
        int max_nwindows
        KmerSlidingWindow* window

    void binary_kmer_alloc_kmers_set(KmerSlidingWindowSet * windows, int max_windows, int max_kmers)
    void binary_kmer_free_kmers_set(KmerSlidingWindowSet **)

    int get_sliding_windows_from_sequence(char*, char*, int, char, short, KmerSlidingWindowSet*, int, int, boolean, int)

###################################################################################################

ctypedef signed char boolean
cdef extern from "element.h":

    ctypedef enum EdgeArrayType:
        pass

    ctypedef enum Element:
        pass

    Key element_get_key(BinaryKmer*, short kmer_size, Key preallocated_key)

###################################################################################################

cdef extern from "hash_table.h":

    ctypedef struct HashTable:
        short kmer_size
        long long number_buckets
        int bucket_size
        Element * table
        short * next_element
        long long * collisions
        long long unique_kmers
        int max_rehash_tries

    HashTable * hash_table_new(int number_bits, int bucket_size, int max_rehash_tries, short kmer_size)
    void hash_table_free(HashTable * * hash_table)
    Element* hash_table_find_or_insert(Key key, boolean * found, HashTable * hash_table)

###################################################################################################

cdef extern from "seq.h":

    ctypedef struct Sequence:
        pass

    void alloc_sequence(Sequence * seq, int max_read_length, int max_name_length)
    void free_sequence(Sequence ** )

###################################################################################################

cdef list assembleReadsAndDetectVariants(char* chrom, int start, int end, cAlignedRead** readStart, cAlignedRead** readEnd, int minQual, int maxReadLength):
    """
    This is the routine which Platypus calls. It creates a dbGraph struct, and sets up a few
    parameters, and then calls the loadAllReadsIntoGraph function to populate the graph. It then
    calls the getVarsFromGraph function to find variant alleles in the graph. The variant alleles
    are returned to Platypus for further processing (mapping, calling etc).
    """
    cdef int bucketSize = 0 #?
    cdef int nBits = 0 #?
    cdef int maxRehashTries = 0 #?
    cdef int kmerSize = 0 #?
    cdef HashTable* graph = hash_table_new(nBits, bucketSize, maxRehashTries, kmerSize)

    #loadAllReadsIntoGraph(graph, readStart, readEnd, minQual, maxReadLength, kmerSize)
    hash_table_free(&graph)

    cdef list cortexVars = []
    cdef Variant v = Variant(chrom, start, "", "", 0, 0, 0, 0, 0)

    return sorted(cortexVars)

###################################################################################################

cdef loadAllReadsIntoGraph(HashTable* graph, cAlignedRead** readStart, cAlignedRead** readEnd, int minQual, int maxReadLength, int kmerSize):
    """
    Take an empty dbGraph struct, and load a set of reads into it. All memory needed for populating
    the graph (except for the graph itself) is allocated and freed here.
    """
    cdef Sequence* seqStruct = <Sequence*>(malloc(sizeof(Sequence)))
    cdef KmerSlidingWindowSet* windows = <KmerSlidingWindowSet*> (malloc(sizeof(KmerSlidingWindowSet)))
    cdef int maxWindows = maxReadLength / (kmerSize+1)
    cdef int maxKmers = maxReadLength - kmerSize + 1

    if seqStruct == NULL or windows == NULL:
        raise StandardError, "Out of memory trying to allocate Sequence or window set"

    alloc_sequence(seqStruct, maxReadLength, LINE_MAX)
    binary_kmer_alloc_kmers_set(windows, maxWindows, maxKmers);

    while readStart != readEnd:
        loadReadIntoGraph(minQual, maxReadLength, graph, readStart[0].seq, readStart[0].qual, seqStruct, windows, kmerSize, maxWindows, maxKmers)

    free_sequence(&seqStruct)
    binary_kmer_free_kmers_set(&windows)

###################################################################################################

cdef loadReadIntoGraph(char minQual, maxReadLength, HashTable* graph, char* seq, char* qual, Sequence* seqStruct, KmerSlidingWindowSet* windows, int kmerSize, int maxWindows, int maxKmers):
    """
    Load a single read into the cortex graph/hash table.
    """
    cdef int index = 0 # TODO: What should this be?
    cdef int entryLength = 0 # TODO: What should this be?
    cdef int markReadStarts = False # TODO: What should this be?
    cdef int fullEntry = True
    cdef int prevFullEntry = True
    cdef int breakHomopolymers = False # TODO: What should this be?
    cdef int homopolymerCutoff = 100 # TODO: What should this be?
    cdef int type = 0 # This can only be 0 (is an enum of EdgeArrayType)

    cdef int nKmers = get_sliding_windows_from_sequence(seq, qual, entryLength, minQual, kmerSize, windows, maxWindows, maxKmers, breakHomopolymers, homopolymerCutoff)
    loadKmersIntoGraph(windows, &prevFullEntry, &fullEntry, markReadStarts, graph, type, index)

###################################################################################################

# A re-write of load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop from "file_reader.c"
cdef loadKmersIntoGraph(KmerSlidingWindowSet* windows, int* prevFullEnt, int* fullEnt, int markReadStarts, HashTable* graph, int type, int index):
    """
    """
    cdef long long totBasesLoaded = 0
    cdef int i = 0
    cdef int j = 0

    cdef Element* currentNode = NULL
    cdef Element* previousNode = NULL
    cdef int currentOrientation = 0 # 0 is forward
    cdef int previousOrientation = 0
    cdef Key tmpKmer
    cdef Key tmpKmer2
    cdef boolean found = False
    cdef KmerSlidingWindow* currentWindow = NULL
    cdef long long lenThisWindow = 0

    # Loop over windows
    for i from 0 <= i < windows[0].nwindows:
        currentWindow = &(windows[0].window[i])
        lenThisWindow = <long long> (currentWindow[0].nkmers + graph[0].kmer_size - 1)

        # Loop over kmers
        for j from 0 <= j < currentWindow[0].nkmers:

            found = False
            tmpKmer2 = element_get_key(&(currentWindow[0].kmer[j]), graph[0].kmer_size, tmpKmer)
            currentNode = hash_table_find_or_insert(tmpKmer2, &found, graph)

            if currentNode == NULL:
                raise StandardError, "file_reader: problem - current kmer not found"
	#  
    #        # Otherwise is the same old last entry
	#        if not (i == 0 && j == 0 && prevFullEnt[0] == false && currentNode == previousNode):
	#            db_node_update_coverage(currentNode, type, index, 1)
	#  
	#        currentOrientation = db_node_get_orientation(&(currentWindow->kmer[j]), currentNode, graph->kmer_size)
	#  
	#        if markReadStarts:

	#            if i == 0 && j == 0 && currentOrientation == forward:
	#	            db_node_set_read_start_status(currentNode, forward)

	#            elif i == 0 && j == 0 && currentOrientation == reverse:
	#	            db_node_set_read_start_status(currentNode, reverse)
	#  
	#        if j > 0:
	#            if previous_node == NULL:
	#                raise StandardError, "file_reader: problem - prev kmer not found"
    #            # This is the point at which the new element/node gets associated with the specific person
	#            else:
	#                db_node_add_edge(previous_node,current_node,previous_orientation,current_orientation, db_graph->kmer_size, type, index)

	#        previous_node = current_node
	#        previous_orientation = current_orientation

###################################################################################################
