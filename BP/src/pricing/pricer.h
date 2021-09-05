#ifndef PRICER_H
#define PRICER_H

#include <vector>
#include <math.h>
#include <cmath>
#include <math.h>       /* sqrt */
#include <numeric>      // std::iota
#include <time.h>
#include <sys/time.h>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <chrono>
#include "../cg.h"
#include "mwisp_inst.h"

namespace GCP{
using namespace std;
class Pricer {

public:
    const double THRESHOLD = -1e-12;
    double cutoff;
    vector<double> dual_values;
    const double INF = 1e8;
    double best_rc = 0.;
    double best_obj=0.;
    long nb_node;
    long nb_edge;
    vector<vector<int>> neg_rc_cols;
    vector<double> neg_rc_vals;
    std::vector<long> sorted_indices;

    virtual ~Pricer(){};
    virtual void run(){std::cout << "error, in pricer!!!\n\n";};

    vector<long> sort_indexes_inc(const vector<double> &v);
    long current_time_for_seeding();
    vector<long> sort_indexes_dec(const vector<double> &v);
    double get_wall_time();
    virtual void add_columns(SCIP_PRICERDATA* pricerdata, bool is_root);
};

}

#endif