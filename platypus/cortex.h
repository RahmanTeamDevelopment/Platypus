#include <binary_kmer.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <dB_graph.h>
#include <seq.h>
#include <file_reader.h>
#include <global.h>
#include <string.h>
#include <dB_graph_supernode.h>
#include <dB_graph_population.h>
#include <file_format.h>
#include <unistd.h>
#include <samtoolsWrapper.c>

void getVarsFromGraph(int maxLength, dBGraph* graph, void (*action_branches)(dBNode*), void (*action_flanks)(dBNode*), Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*));
void get_vars(dBNode* node);
void get_vars_with_orientation(dBNode* node, Orientation orientation);
