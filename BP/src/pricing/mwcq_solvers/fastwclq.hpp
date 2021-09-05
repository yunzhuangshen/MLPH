#ifndef FASTWCLQ_H
#define FASTWCLQ_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <map>
#include <algorithm>
#include <set>
#include <string.h>

namespace FASTWCLQ_SOLVER {
using namespace std;

struct Remaining_vertex {
	vector<long> vertex;
	vector<vector<long>::size_type> index;

	vector<long>::iterator begin() {
		return vertex.begin();
	}
	vector<long>::iterator end() {
		return vertex.end();
	}
	void init(vector<long>::size_type vertex_size) {
		vertex.reserve(vertex_size);
		index.resize(vertex_size);
		for (vector<long>::size_type i = 1; i < vertex_size; ++i) {
			vertex.push_back(i);
			index[i] = i - 1;
		}
	}
	void remove(long v) {
		index[*vertex.rbegin()] = index[v];
		vertex[index[v]] = *vertex.rbegin();
		vertex.pop_back();
	}

	vector<long>::size_type size() {
		return vertex.size();
	}

	bool empty() {
		return vertex.empty();
	}
};

class FASTWCLQ{

const float       MY_RAND_MAX_FLOAT = 10000000.0;
const long long   	  MY_RAND_MAX_int =   10000000;
const float 	  BASIC_SCALE = 0.0000001; //1.0f/MY_RAND_MAX_FLOAT;
vector<vector<long>> adjacency_list;
vector<long> vertex_neighbor_weight;
vector<long> vertex_weight;
//reduction
vector<long> hit_in_common_neighbor;
vector<long> vertex_to_removed;
vector<long> working_vertex;
vector<long> next_working_vertex;
vector<bool> is_pending;
vector<vector<long>::size_type> index_in_working_vertex;
//solutions
vector<long> solution;
vector<long> best_solution;
long long best_solution_weight=0;
long long solution_weight;
long tries;
double best_solution_time;
long best_solution_try;
//input parameter
long size_threshold;
long long t;
Remaining_vertex remaining_vertex;
vector<long> start_vertices;
long untest_pointer;
vector<vector<long long>> adjacency_cand_neighbor_weight;
vector<bool> is_computed; 
vector<long> candidates;
vector<long long> cand_neighbor_weight;
vector<bool> is_in_candidates;
vector<bool> is_addv_neighbor; // indicates whether a candidate vertex is adjacent to the add_v
//vector<long> remove_cand_vertices;//in each step
//vector<bool> is_removed;
long start_bms_count=1;
long min_bms_count;
long max_bms_count;
long real_bms_count;
bool is_new_graph = true;
long simp_count=0;
bool better_since_simp=true;


void build(string file_name);
void print_vec(vector<long> & container) ;
long upper_bound(long v);
void simplify_iterative();
void simplify();
template<typename T>
bool is_neighbor(T v1, T v2);
void output_graph_size();
void output_best_solution();
void update_best_solution();
bool compare_vertex (long v1, long v2);
void init();
long long construct();
vector<long>::size_type binary_search(const vector<long> &array, long value);
bool verify_simple(string file_name);
bool verify (string file_name);


public:
	FASTWCLQ(){};
    std::vector<double> sol_objs;
    std::vector<double> sol_times;
    std::vector<std::vector<int>> sols;
    long long fastwclq(long ntries, long seed, long long NB_Node, long long NB_EDGE, double cutoff, long long** AdjacentList, long long* Node_Degree, long long* Node_Weight, long long* Node_Bound);

};
}
#endif