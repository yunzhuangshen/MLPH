#ifndef TSM_H
#define TSM_H


#include <vector>
#include <iostream>


namespace TSM {
    

#define WORD_LENGTH 100
#define TRUE 1
#define FALSE 0
#define NONE -1
#define DELIMITER 0
#define PASSIVE 0
#define ACTIVE 1
#define UNASSIGNED 2
#define P_TRUE 1
#define P_FALSE 0
#define NULL_REASON -1
#define NO_REASON -3
#define CONFLICT -1978
#define MAX_NODE 1000000
#define max_expand_depth 100000
#define STACK_LENGTH (MAX_NODE*2)
#define pop(stack) stack[--stack ## _fill_pointer]
#define push(item, stack) stack[stack ## _fill_pointer++] = item
#define ptr(stack) stack ## _fill_pointer
#define is_neibor(i,j) matrice[i][j]

#define CUR_CLQ_SIZE Clique_Stack_fill_pointer
#define CURSOR Cursor_Stack[Cursor_Stack_fill_pointer-1]
#define MIN(a,b) a<=b?a:b
#define BIT_MAP_SIZE 4097
#define Node_Reason Node_Degree

#define SET_EDGE(row,col) ((*(Adj_Matrix + (row)* MATRIX_ROW_WIDTH + ((col) >> 3))) |= (1 << ((col) & 7)))
#define GET_EDGE(row,col) ((*(Adj_Matrix + (row)* MATRIX_ROW_WIDTH + ((col) >> 3))) & (1 << ((col) & 7)))

#define iMatrix(i) (Adj_Matrix+(i)*MATRIX_ROW_WIDTH)
#define Matrix(i,j) ((*((i) + ((j) >> 3))) & (1 << ((j) & 7)))
#define New_Name Node_Degree
#define assign_node(node, value, reason) {\
	Node_Value[node] = value;\
	Node_Reason[node] = reason;\
	push(node, FIXED_NODE_STACK);\
}
struct iSET_State {
	char satisfied;
	char used;
	char involved;
	char active;
	int size;
	int topk;
	long long weight;  //int unassigned;
	long long t_weight;  //int unassigned;
	long long *nodes;
};

class TSM{
int * Adj_List;

int FORMAT = 1, NB_NODE, NB_NODE_O, NB_EDGE, NB_EDGE_O, MAX_CLQ_SIZE, INIT_CLQ_SIZE, INIT_ORDERING, NB_BACK_CLIQUE, MATRIX_ROW_WIDTH, MATRIX_SIZE = 0, MAX_SUBGRAPH_SIZE, K_CORE_G = 0;

long long MAX_CLQ_WEIGHT, CUR_CLQ_WEIGHT, OPT_CLQ_WEIGHT, MAX_ISET_WEIGHT, CUR_NODE_WEIGHT, UPPER_WEIGHT_BOUND;
long long INIT_CLQ_WEIGHT = 0;
long long Max_Degree = 0;

long long* Node_Degree;
long long* Top_Weight;
long long* Node_Weight;
char* Node_Value;
int* rankVar;
long long* Candidate_Stack; //MAX_NODE * 2
long long* Vertex_UB; //MAX_NODE * 2

int **Node_Neibors;

int Candidate_Stack_fill_pointer = 0;
int Clique_Stack_fill_pointer;
int *Clique_Stack, *MaxCLQ_Stack;
int Cursor_Stack[max_expand_depth];
int Cursor_Stack_fill_pointer = 0;
int Weight_Mod = 200;
unsigned char * Adj_Matrix;

int iSET_COUNT = 0;
long long iSET_TOTAL_WEIGHT = 0;
int *iSET_Size;
long long *iSET_Weight;
char *iSET_Tested;
long long **iSET;
int NEW_NODE_IDX = 0, MAX_OLD_NODE = 0, MAX_ISET_COUNT = 0;
int RESERVED_LENGTH = 100;
int LARGE_WEIGHT = FALSE;

int threshold;
char* File_Name;
const char* resPath;
const char* dataPath;

/* the structures for maxsat reasoning*/
struct iSET_State *IS;
long long *iNode_TAIL;
long long *REASON_STACK;
int REASON_STACK_fill_pointer = 0;
int *UNIT_STACK;
int UNIT_STACK_fill_pointer = 0;
//int *TOP_UNIT_STACK;
//int TOP_UNIT_STACK_fill_pointer = 0;
int *NEW_UNIT_STACK;
int NEW_UNIT_STACK_fill_pointer = 0;
int *FIXED_NODES_STACK;
int FIXED_NODES_STACK_fill_pointer = 0;
int *SATISFIED_iSETS_STACK;
int SATISFIED_iSETS_STACK_fill_pointer = 0;
int *TOPK_REDUCED_STACK;
int TOPK_REDUCED_STACK_fill_pointer = 0;

int Rollback_Point;
int Branching_Point;
int *Old_Name;
int *Second_Name;
int NB_CANDIDATE = 0, FIRST_INDEX;
int START_MAXSAT_THD = 15;
int BRANCHING_COUNT = 0;
int CUR_MAX_NODE;

int Last_Idx = 0;
int cut_ver = 0, total_cut_ver = 0;
int cut_inc = 0, total_cut_inc = 0;
int cut_iset = 0, total_cut_iset = 0;
int cut_satz = 0, total_cut_satz = 0;
int MATRIX_REBUILDED = 0;
int REBUILD_MATRIX = FALSE;
int MATRIX_LIMITATION = 2048;

int LAST_IN;
int NB_TOPK = 0;
/*statistical data*/
long long N0_B, N0_A, N1_B, G1_B = 0;
double D0_B, D0_A, D1_B = 0;
/*****************/
int * Init_Adj_List;
int BLOCK_COUNT = 0;
int *BLOCK_LIST[100];
double READ_TIME, INIT_TIME, SEARCH_TIME;
int NODE_IN_ISET;
double get_utime();
void free_block();
int is_adjacent(int node1, int node2);
void check_clique_in_result_file(char *input_file);
void allcoate_memory_for_adjacency_list(int nb_node, int nb_edge, int offset);
int read_graph_wclq_format(char *input_file);
int read_graph_node_node(char *input_file, int format);
int build_simple_graph_instance(char *input_file);
void sort_by_degree_degeneracy_ordering();
void sort_by_score_ordering();
long long init_for_degeneracy(long long *segment_counter, long long *where);
long long sort_by_weight_degeneracy_and_core_decomposion();
int addIntoIsetTomitaBis_adj(int node);
int re_number_adj(int node);
long long absorb_by_inserting(int iset_idx, long long iset_weight, int insert_node, long long node_weight, int last_iset);
long long absorb_by_splitting(int iset_idx, long long topk_weight, int insert_node, long long node_weight);
void do_weight_partition(int insert_node, long long node_weight);
int add_vertex_with_weight_partition(int insert_node);
int cut_by_iset_less_vertices();
int cut_by_binary_maxsat_reasoning();
int add_new_iset_for_bnode(int b_node, long long b_weight);
void init_for_maxsat_reasoning();
int retrace_conflicting_isets(int start_iset);
int fix_node(int fix_node, int fix_iset);
int reasoning_with_up();
void reset_context();
void reduce_iset(int iset_idx, long long delta);
void reduce_iset_topk_new(int iset_idx, long long delta, long long t_weight);
void split_conflicting_isets();
long long insert_into_iset(int iset_idx, long long iset_weight, int insert_node, long long node_weight);
long long split_and_insert(int iset_idx, long long t_weight, int insert_node, long long node_weight);
int insert_vertex_with_max3sat(int insert_node);
int cut_by_ordered_maxsat_reasoning();
void build_iset_structures();
void build_node_isets();
void rebuild_matrix(int start);
void store_maximal_weighted_clique();
long long compute_subgraphs_weight(int start);
void store_maximal_weighted_clique2();
long long init_for_me(int start, long long *segment_counter, long long *where);
int compute_subgraph_degree(int start);
int reduce_first_level_subgraphs_by_degree(int start);
int reduce_first_level_subgraphs_by_weight(int start);
int cut_by_inc_ub();
void init_for_search();
void allocate_memory();
void search_maxclique(double cutoff, int using_init_clique);
void check_maxw_clique();
void printMaxClique(char* input_file);
void build_init_matrix();
void reduce_instance_for_weight();
void reduce_instance_for_degree();
void print_version1();
void check_result(char *input_file, char *result_file);
void clear_structures();


public:
TSM();
~TSM();
bool isOptimal;
std::vector<double> sol_objs;
std::vector<double> sol_times;
std::vector<std::vector<int>> sols;
long long tsm(char* File_Name, const char* resPath, int NB_Node, int NB_EDGE, double cutoff, long long** AdjacentList, long long* Node_Degree, long long* Node_Weight, long long* Node_Bound);


};

}
#endif