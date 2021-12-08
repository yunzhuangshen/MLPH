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
#include "pricer.h"


namespace GCP {
  using namespace std;
  class Heur_pricer_misc : public Pricer{
    
  public:
    METHOD_TYPE method_type;
    const vector<vector<int>>& mis_adj_list;
    const vector<vector<bool>>& mis_adj_matrix;
    
    vector<double> objs;
    vector<vector<int>> sols;

    explicit Heur_pricer_misc(METHOD_TYPE _method_type, double _cutoff, MWISP_INST& inst);
    void solve_mwis_greedy();
    void run() override;
  };
}

#endif
