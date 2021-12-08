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

namespace GCP{
using namespace std;
class Pricer {

public:
    int method;
    int column_selection = 0;
    const double THRESHOLD = -0.000001;
    double cutoff;
    int sample_size;
    int niterations;
    int basis_factor;
    int upper_col_limit;

    vector<double> dual_values;
    double heur_best_reduced_cost;
    const double INF = 1e8;

    long nb_node;
    long nb_edge;
    vector<vector<int>> neg_rc_cols;
    vector<double> neg_rc_vals;
    std::vector<long> sorted_indices;
    long num_neg_rc_col=0;;
    double best_rc = 0.;
    double mean_rc = 0.;
    double stdev_rc = 0.;
    double median_rc = 0.;

    //only for mssp-all and mssp-clustering 
    double column_selection_time=0.;
    //only for mssp-all and mssp-dominance selection
    double dominate_selection_time=0;

    vector<double> compute_nrc(vector<vector<int>>& basic_cols);

    virtual ~Pricer(){};
    virtual void run(){std::cout << "error, in pricer!!!\n\n";};

    vector<long> sort_indexes_inc(const vector<double> &v);
    long current_time_for_seeding();
    vector<long> sort_indexes_dec(const vector<double> &v);

    double get_wall_time();

    void compute_statistics();

    void include_new_cols(std::vector<std::vector<int>>& basic_cols, vector<int>& lb_vbasis);
    
    void add_all(std::vector<std::vector<int>>& basic_cols);
    void add_partial(std::vector<std::vector<int>>& basic_cols);
    void replace_existing(vector<vector<int>>& basic_cols, vector<int>& lb_vbasis );

};

}

#endif