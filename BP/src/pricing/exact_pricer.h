#ifndef Exact_pricer_H
#define Exact_pricer_H
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
#include "../BP/cons_storeGraph.h"

namespace GCP {
  using namespace std;
  class Exact_pricer : public Pricer{
    
  public:
    double best_pricing_obj;
    double cutoff;
    int nb_node;
    EXACT_PRICER exact_type;
    const vector<vector<int>>& mis_adj_list;
    const vector<vector<bool>>& mis_adj_matrix;
    MWISP_INST& inst;
    SCIP_PRICERDATA* pricerdata;
    vector<double> objs;
    vector<double> times;
    vector<vector<int>> sols;
    vector<double> reduced_costs;
    bool isOptimal;
    double max_mwis_obj;
    // Builds a Exact_pricer for graph 
    explicit Exact_pricer(EXACT_PRICER _exact_type, double cutoff, MWISP_INST& inst, SCIP_PRICERDATA* pricerdata, double best_pricing_obj);
    void solve_tclique();
    void solve_tsm();
    void solve_gurobi();
    void run() override;
  };
}

#endif
