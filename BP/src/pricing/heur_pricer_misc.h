#ifndef Heur_pricer_misc_H
#define Heur_pricer_misc_H
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <numeric>      // iota
#include <cassert>
#include <set>
#include <gurobi_c++.h>
#include "mwcq_solvers/lscc.hpp"
#include "mwcq_solvers/fastwclq.hpp"
#include "mwcq_solvers/tsm.hpp"
#include "pricer.h"


namespace GCP {
  using namespace std;
  class Heur_pricer_misc : public Pricer{
    
  public:
    METHOD_TYPE method_type;
    const vector<vector<int>>& mis_adj_list;
    const vector<vector<bool>>& mis_adj_matrix;
    long long** AdjacentList;
    long long* Node_Degree;
    long long* Node_Weight;
    long long* Node_Bound;
    
    vector<double> objs;
    vector<vector<int>> sols;

    explicit Heur_pricer_misc(METHOD_TYPE _method_type, double _cutoff, MWISP_INST& inst);
    void solve_mwis_greedy();
    void solve_mwc_tsm();
    void solve_mwc_wlmc();
    void solve_mwc_lscc();
    void solve_mwc_fastwclq();
    void solve_mwis_gurobi();
    void run() override;
    void postprocessing();
    ~Heur_pricer_misc();
  };
}

#endif
