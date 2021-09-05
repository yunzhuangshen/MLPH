#ifndef Random_col_H
#define Random_col_H
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <numeric>      // iota
#include <cassert>
#include <set>
#include "pricer.h"
#include "../BP/cons_storeGraph.h"
#include <random>

namespace GCP {
  using namespace std;
  class Random_col{
    
  public:
    int nb_node;
    int seed;
    int nb_mwis;
    const vector<vector<bool>>& mis_adj_matrix;
    explicit Random_col(MWISP_INST& inst, int seed, int nb_mwis);
    void add_random_cols(SCIP* scip, SCIP_CONS** constraints);
  };
}

#endif
