#include "cortex.h"

//-------------------------------------------------------------------------------------------------
//
// TODO: 1) Initialise graph
//       2) Add reads to graph (set max chunk length to large value e.g. 10,000 for reference reading)
//       3) Clean graph
//       4) Add reference to graph
//       5) Find bubbles in graph
//       6) Call variants
//
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
//routine to DETECT/DISCOVER variants directly from the graph - reference-free (unless you have put the reference in the graph!)
// "condition" argument is a condition which you apply to the flanks and branches to decide whether to call.
// e.g. this might be some constraint on the coverage of the branches, or one might have a condition that one branch
//      was in one colour  and the other in a different colour, or maybe that both branches are in the same colour
// last argument get_colour specifies some combination of colours, defining the graph within which we look for bubbles.
// most obvious choices are: colour/edge with index (say)2, or union of all edges, or union of ll except that which is the reference genome
//-------------------------------------------------------------------------------------------------

void getVarsFromGraph(int maxLength, dBGraph* graph, void (*action_branches)(dBNode*), void (*action_flanks)(dBNode*), Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*))
{
   // int count_vars = 0; 
   // int flanking_length = 1000; 

   // // Allocate arrays
   // dBNode** path_nodes1 = (dBNode**) malloc(sizeof(dBNode*)*(max_length+1));
   // dBNode** path_nodes2 = (dBNode**) malloc(sizeof(dBNode*)*(max_length+1));
   // Orientation* path_orientations1 = (Orientation*) malloc(sizeof(Orientation)*(max_length+1));
   // Orientation* path_orientations2 = (Orientation*) malloc(sizeof(Orientation)*(max_length+1));
   // Nucleotide* path_labels1 = (Nucleotide*) malloc(sizeof(Nucleotide)*(max_length+1) );
   // Nucleotide* path_labels2 = (Nucleotide*) malloc(sizeof(Nucleotide)*(max_length+1) );
   // char* seq1 = (char*) malloc(sizeof(char)*(max_length+1));
   // char* seq2 = (char*) malloc(sizeof(char)*(max_length+1));

   // if ( (path_nodes1==NULL) || (path_nodes2==NULL) || (path_orientations1==NULL) || (path_orientations2==NULL) 
   //         || (path_labels1==NULL) || (path_labels2==NULL) || (seq1==NULL) || (seq2==NULL) )
   // {
   //     printf("Could not allocate arrays in db_graph_detect_vars. Out of memory, or you asked for unreasonably big max branch size: %d\n", max_length);
   //     exit(1);
   // }

   // hash_table_traverse(&get_vars,db_graph); 

   // // Cleanup
   // free(path_nodes1);
   // free(path_nodes2);
   // free(path_orientations1);
   // free(path_orientations2);
   // free(path_labels1);
   // free(path_labels2);
   // free(seq1);
   // free(seq2);
}

//-------------------------------------------------------------------------------------------------

void get_vars(dBNode* node)
{
   // if (db_node_check_status_none(node))
   // {     
   //     action_flanks(node);
   //     get_vars_with_orientation(node,forward);
   //     get_vars_with_orientation(node,reverse);
   // }
}

//-------------------------------------------------------------------------------------------------

void get_vars_with_orientation(dBNode* node, Orientation orientation)
{
    //int length1, length2;
    //double avg_coverage1;
    //double avg_coverage2;
    //int min_coverage1,max_coverage1;
    //int min_coverage2,max_coverage2;

    //dBNode* current_node = node;
    //VariantBranchesAndFlanks var; //will reuse this as we traverse the graph

    //do
    //{
    //    // The idea is that db_graph_detect_bubble_in_subgraph_defined_by_func_of_colours will mark anything is sees with action_flanks
    //    // However if it does find a bubble, and if you want it to (penultimate argument=true)  then the branches are  marked with action_branches

    //    if (db_graph_detect_bubble_in_subgraph_defined_by_func_of_colours(current_node,orientation,max_length,action_flanks,
    //                &length1,path_nodes1,path_orientations1,path_labels1,
    //                seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
    //                &length2,path_nodes2,path_orientations2,path_labels2,
    //                seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
    //                db_graph, get_colour, get_covg,
    //                true, action_branches))
    //    {

    //        int length_flank5p = 0;	
    //        int length_flank3p = 0;
    //        dBNode * nodes5p[flanking_length];
    //        dBNode * nodes3p[flanking_length];
    //        Orientation orientations5p[flanking_length];
    //        Orientation orientations3p[flanking_length];
    //        Nucleotide labels_flank5p[flanking_length];
    //        Nucleotide labels_flank3p[flanking_length];
    //        char seq5p[flanking_length+1];
    //        char seq3p[flanking_length+1];
    //        boolean is_cycle5p, is_cycle3p;

    //        double avg_coverage5p;
    //        int min5p,max5p;
    //        double avg_coverage3p;

    //        int min3p,max3p;

    //        //compute 5' flanking region       	
    //        int length_flank5p_reverse = db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(current_node,opposite_orientation(orientation),
    //                flanking_length,action_flanks,
    //                nodes5p,orientations5p,labels_flank5p,
    //                seq5p,&avg_coverage5p,&min5p,&max5p,
    //                &is_cycle5p,db_graph, get_colour, get_covg);

    //        if (length_flank5p_reverse>0){
    //            Nucleotide label;
    //            Orientation next_orientation;

    //            dBNode * lst_node = db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(nodes5p[length_flank5p_reverse-1],orientations5p[length_flank5p_reverse-1],
    //                    &next_orientation, labels_flank5p[length_flank5p_reverse-1],&label,
    //                    db_graph, get_colour);

    //            length_flank5p = db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(nodes5p[length_flank5p_reverse],
    //                    opposite_orientation(orientations5p[length_flank5p_reverse]),
    //                    flanking_length,label,
    //                    action_flanks,
    //                    nodes5p,orientations5p,labels_flank5p,
    //                    seq5p,&avg_coverage5p,&min5p,&max5p,
    //                    &is_cycle5p,db_graph, get_colour, get_covg);
    //        }
    //        else{
    //            length_flank5p = 0;
    //        }

    //        //compute 3' flanking region
    //        length_flank3p = db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(path_nodes2[length2],path_orientations2[length2],
    //                flanking_length, action_flanks,
    //                nodes3p,orientations3p,labels_flank3p,
    //                seq3p,&avg_coverage3p,&min3p,&max3p,
    //                &is_cycle3p,db_graph, get_colour, get_covg);

    //        char name[100];

    //        //warning - array of 5prime nodes, oprientations is in reverse order to what you would expect - it is never used in what follows
    //        set_variant_branches_and_flanks(&var, nodes5p, orientations5p, length_flank5p, path_nodes1, path_orientations1, length1, 
    //                path_nodes2, path_orientations2, length2, nodes3p, orientations3p, length_flank3p, unknown);

    //        if (condition(&var)==true) 
    //        {
    //            count_vars++;

    //            //print flank5p - 
    //            sprintf(name,"var_%i_5p_flank",count_vars);

    //            print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(fout,name,length_flank5p,avg_coverage5p,min5p,max5p,
    //                    nodes5p[0],orientations5p[0],			
    //                    nodes5p[length_flank5p],orientations5p[length_flank5p],				
    //                    seq5p,
    //                    db_graph->kmer_size,true, get_colour, get_covg);	

    //            //print branches
    //            sprintf(name,"branch_%i_1",count_vars);
    //            print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(fout,name,length1,
    //                    avg_coverage1,min_coverage1,max_coverage1,
    //                    path_nodes1[0],path_orientations1[0],path_nodes1[length1],path_orientations1[length1],
    //                    seq1,
    //                    db_graph->kmer_size,false, get_colour, get_covg);

    //            sprintf(name,"branch_%i_2",count_vars);
    //            print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(fout,name,length2,
    //                    avg_coverage2,min_coverage2,max_coverage2,
    //                    path_nodes2[0],path_orientations2[0],path_nodes2[length2],path_orientations2[length2],
    //                    seq2,
    //                    db_graph->kmer_size,false, get_colour, get_covg);

    //            //print flank3p
    //            sprintf(name,"var_%i_3p_flank",count_vars);
    //            print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(fout,name,length_flank3p,avg_coverage3p,min3p,max3p,
    //                    nodes3p[0],orientations3p[0],nodes3p[length_flank3p],orientations3p[length_flank3p],
    //                    seq3p,
    //                    db_graph->kmer_size,false, get_colour, get_covg);
    //        }

    //        action_branches(path_nodes2[length2]);
    //        current_node = path_nodes2[length2];
    //        orientation = path_orientations2[length2];
    //    }
    //} while (current_node != node //to avoid cycles
    //        && db_node_check_status_none(current_node));
}

//-------------------------------------------------------------------------------------------------
