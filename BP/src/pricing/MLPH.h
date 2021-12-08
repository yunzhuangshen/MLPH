#ifndef MLPH_H
#define MLPH_H

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>
#include <set>
#include <cstring>
#include <string>
#include <iomanip>
#include "gurobi_c++.h"
#include "pricer.h"
using namespace std;

namespace GCP{
    class MLPH: public Pricer{
        std::vector<float> predicted_value;
        vector<vector<int>> mis_set;
        const vector<vector<bool>> adj_matrix;
        const vector<double> degree_norm;
        std::vector<double> objs;
        vector<double> dists;

        public:
            METHOD_TYPE variant;
            int random_sample_size;
            int sample_size;
            int niterations;
            double start_time;
            float min_cbm, max_cbm;
            float min_rbm, max_rbm;
            double max_dual;
            std::vector<float> ranking_scores;
            std::vector<float> corr_xy;
            std::vector<double> bound_norm;
            double b0;
            double b1;
            double c0;
            double c1; 
            double c2; 
            double c3; 
            double c4;
            double b;
            std::set<std::string> identites;
            MLPH(bool is_root, METHOD_TYPE variant, double cutoff, MWISP_INST& inst, double sample_factor);
            void run();
            void random_sampling();
            void make_prediction();            
            void _run();           
            void compute_correlation_based_measure();
            void compute_ranking_based_measure();    
            void compute_bound();
    };
}

#endif
