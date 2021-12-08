#ifndef CG_H
#define CG_H
#include "instance.h"

extern "C" {
#include "svm_predict_model.h"
#include "linear_svm_predict_model.h"
}

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <numeric>      // std::iota
#include <vector>
#include <cstring>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <iomanip>
#include "gurobi_c++.h"

using namespace std;

namespace GCP{
    class CG{
        const Instance& g;
        const string test_data_dir = "../test_data/";
        const string training_data_dir = "../train_data/";
        const string training_model_name = "train_model";
        const string test_data_name = "_test_data_";
        const string output_file_name = "_predicted_value_";
        const int mis_factor = 10;
        double min_reduced_cost_exact;

        vector<vector<int>> mis_set;
        vector<double> dual_value;
        vector<vector<int>> adj_list;
        vector<vector<bool>> adj_matrix;
        vector<int> optimal_mis;

        void initializing_parameters();
    public:
        int seed;
        double time_computing_lagrangian_bound=0;
        double cutoff;
        double pricer_cutoff;
        int thread_limit=1;
        int num_heur_runs_success=0;
        double time_duration_master=0;
        double time_duration_pricing_exact=0;
        double time_duration_pricing_heur=0;
        bool lp_optimal=false;
        vector<int> lp_vbasis;
        vector<int> lp_cbasis;
        double lp_bound = 0.0;
        int cg_iters = 0;
        int num_mis = 0;

        explicit CG(const Instance& g, double cutoff, double pricer_cutoff, double _thread_limit, int _seed);
        void initializing_mis();
        void solve_restricted_master_problem(bool warm_start);
        bool solve_mwis_gurobi(double cutoff, double& min_rc);
        bool solve_mwis_tsm(double cutoff, double& min_rc);

        void collect_training_data(vector<vector<double>>& obj_coef, vector<vector<bool>>& solution);
        double optimize_LM(int method, double b0, double b1);
        void test(int method, int column_selection, std::ofstream* output_file_cg_stats=nullptr);        
    };
}

#endif
