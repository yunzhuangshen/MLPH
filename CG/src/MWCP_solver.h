#ifndef MWCP_solver_H
#define MWCP_solver_H
#include <iostream>
#include <vector>
#include "mwcq_solvers/lscc.hpp"
#include "mwcq_solvers/fastwclq.hpp"
#include "mwcq_solvers/tsm.hpp"
#include "pricer.h"

namespace GCP {
  using namespace std;
  class MWCP_solver : public Pricer{
    
  public:
    double exact_rc;
    const vector<vector<int>>& mis_adj_list;
    const vector<vector<bool>>& mis_adj_matrix;
    std::vector<long long> nodes_list;
    long long** AdjacentList;
    long long* Node_Degree;
    long long* Node_Weight;
    long long* Node_Bound;
    
    vector<double> objs;
    vector<double> times;
    vector<vector<int>> sols;
      vector<double> reduced_costs;

    //only for tsm
    vector<int> optimal_mis;
    bool isOptimal=false;
    double best_obj;
    // Builds a MWCP_solver for graph 
    explicit MWCP_solver(int _method, double cutoff, const vector<double>& _dual_values, 
        const vector<vector<int>>& _mis_adj_list, const vector<vector<bool>>& mis_adj_matrix, 
        long long nb_edges, long long _upper_col_limit);
    void solve_mwis_greedy();
    void solve_mwc_tsm();
    void solve_mwc_wlmc();
    void solve_mwc_lscc();
    void solve_mwc_fastwclq();
    void solve_mwis_gurobi();
    void run() override;
    void postprocessing();
    ~MWCP_solver();
  };
}

#endif
