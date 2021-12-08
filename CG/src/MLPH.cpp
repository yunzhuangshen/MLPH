#include "MLPH.h"
#include <limits.h>
#include <sstream>
#include <cassert>
#include <omp.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
namespace GCP {


    MLPH::MLPH(int _method, double b0, double b1, int _n, int _sample_size, const vector<double>& _degree_norm, 
                const vector<vector<bool>>& _adj_matrix, const vector<double>& _dual_values,
                int _upper_col_limit): degree_norm{_degree_norm}, adj_matrix{_adj_matrix}, b0{b0}, b1{b1}{
        
        dual_values = _dual_values;
        method = _method;
        nb_node=_n;
        sample_size = _sample_size;
        upper_col_limit = _upper_col_limit;
        niterations = 1;
        heur_best_reduced_cost = 1;
        assert(method==6);

        c0=1.6557295690309566;
        c1=-1.061904300496613;
        c2=-4.631979340281991;
        c3=-1.5341779416628611;
        c4=5.4064285368308;

        // read svm.param if exists
        std::ifstream svm_param_file("../../svm.param");
        if (svm_param_file.good()){
            string line;
            getline(svm_param_file, line);
            c0 = std::stod(line);
            getline(svm_param_file, line);
            c1 = std::stod(line);
            getline(svm_param_file, line);
            c2 = std::stod(line);
            getline(svm_param_file, line);
            c3 = std::stod(line);
            getline(svm_param_file, line);
            c4 = std::stod(line);
            getline(svm_param_file, line);
            b = std::stod(line);
            svm_param_file.close();
        }

        cutoff = 1e10;
        max_dual = 0.;
        for (auto i = 0; i < nb_node; i++){
            // cout << i << " " << dual_values[i] << "\n";
            if (dual_values[i] > max_dual)
                max_dual = dual_values[i];
        }
        compute_bound();
        identites = std::set<std::string>();
    }

    MLPH::MLPH(int _method, double _cutoff, int _n, int _sample_size, const vector<double>& _degree_norm, 
                const vector<vector<bool>>& _adj_matrix, const vector<double>& _dual_values,
                int _upper_col_limit): 
        degree_norm{_degree_norm}, adj_matrix{_adj_matrix}{
            
        dual_values = _dual_values;
        method = _method;
        assert(method==6);
        nb_node = _n;
        cutoff = _cutoff;
        sample_size = _sample_size;
        upper_col_limit = _upper_col_limit;
        niterations = 50;
        heur_best_reduced_cost = 1;

        c0=1.6557295690309566;
        c1=-1.061904300496613;
        c2=-4.631979340281991;
        c3=-1.5341779416628611;
        c4=5.4064285368308;
        //1-1.0/(1+e^(9.775*x-12.5564))
        b0 = -9.775;
        b1 = -12.5564;

        // read svm.param if exists
        std::ifstream svm_param_file("../../svm.param");
        if (svm_param_file.good()){
            string line;
            getline(svm_param_file, line);
            c0 = std::stod(line);
            getline(svm_param_file, line);
            c1 = std::stod(line);
            getline(svm_param_file, line);
            c2 = std::stod(line);
            getline(svm_param_file, line);
            c3 = std::stod(line);
            getline(svm_param_file, line);
            c4 = std::stod(line);
            getline(svm_param_file, line);
            b = std::stod(line);
            svm_param_file.close();
            // std::cout << c0 << " " << c1 << " " << c2 << " " << c3 << " " << c4 << " " << b << "\n"; 
        }
        // read lm.param if exists
        std::ifstream lm_param_file("../../lm.param");
        if (lm_param_file.good()){
            string line;
            getline(lm_param_file, line);
            b0 = std::stod(line);
            getline(lm_param_file, line);
            b1 = std::stod(line);
            lm_param_file.close();
            // std::cout << b0 << " " << b1 << "\n";
        }

        max_dual = 0.;
        for (auto i = 0; i < nb_node; i++){
            // cout << i << " " << dual_values[i] << "\n";
            if (dual_values[i] > max_dual)
                max_dual = dual_values[i];
        }
        compute_bound();
        identites = std::set<std::string>();
    }

    void MLPH::random_sampling() {
        objs = std::vector<double> (sample_size);
        mis_set = std::vector<std::vector<int>> (sample_size);

        long time_seed = current_time_for_seeding();        
        #pragma omp parallel
        {
            int threadnum = omp_get_thread_num();
            mt19937 mt(time_seed + threadnum*1e5);
            uniform_int_distribution<int> dist(0,RAND_MAX);

            int v, idx, num;
            vector<int> candidates(nb_node);
            int nb_candidates;
            #pragma omp for
            for(int i = 0; i < sample_size; ++i) {
                objs[i]=0;
                nb_candidates = nb_node;
                for (int j = 0; j < nb_candidates; ++j){
                    candidates[j] = j;
                }
                while (nb_candidates > 0){
                    if (nb_candidates == nb_node){
                        idx = i % nb_node;
                    } else{
                        idx = dist(mt) % nb_candidates;
                    }
                    v = candidates[idx];
                    mis_set[i].push_back(v);
                    objs[i] += dual_values[v];
                    
                    num = 0;
                    for (int j = 0; j < nb_candidates; ++j){
                        if (adj_matrix[v][candidates[j]] == 0 && j != idx){
                            candidates[num] = candidates[j];
                            num++;
                        }
                    }
                    nb_candidates = num;
                }

                std::sort (mis_set[i].begin(), mis_set[i].end());  
                std::stringstream ss;
                for (auto k = 0; k < mis_set[i].size(); k++){
                    ss << mis_set[i][k] << " ";
                }
                std::string identity = ss.str();  
                #pragma omp critical
                {
                    if (identites.find(identity) == identites.end()){
                        if (1 - objs[i] < THRESHOLD){
                            neg_rc_cols.push_back(mis_set[i]);
                            neg_rc_vals.push_back(1-objs[i]);
                        }
                        identites.insert(identity);
                    }
                }
            }
        }
    }

    void MLPH::compute_bound(){
        bound_norm = std::vector<double> (nb_node, 0);
        double max_bound=0;
        double min_bound=1e8;
        for (auto i = 0; i < nb_node; i++){
            for (auto j = 0; j < nb_node; j++){
                if (i!=j && adj_matrix[i][j] == 0)
                    bound_norm[i] += dual_values[j]; 
            }
            if (bound_norm[i] > max_bound)
                max_bound = bound_norm[i];
            if (bound_norm[i] < min_bound)
                min_bound = bound_norm[i];
        }

        double bound_delta = max_bound - min_bound;
        if (bound_delta != 0){
            for (auto i = 0; i < nb_node; i++)
                bound_norm[i] /= max_bound;
        }
    }
    void MLPH::compute_ranking_based_measure() {

        std::vector<int> sort_idx(sample_size);
        std::iota(sort_idx.begin(), sort_idx.end(), 0);
        std::vector<double> objs_copy(objs);
        std::sort(sort_idx.begin(), sort_idx.end(), [&objs_copy](int i1, int i2) {return objs_copy[i1] > objs_copy[i2];});
        std::vector<int> rank = std::vector<int>(sample_size);
        for (auto i = 0; i < sample_size; i++){
            rank[sort_idx[i]] = i;
        }
        
        std::vector<int> sampling_ct(nb_node, 0);
        ranking_scores = std::vector<float>(nb_node, 0.0);
        for (auto i = 0; i < sample_size; ++i){
            for (auto j = 0; j < mis_set[i].size(); ++j){
                sampling_ct[mis_set[i][j]]++;
                ranking_scores[mis_set[i][j]] += 1.0/(rank[i]+1);
            }
        }

        min_rbm=1.; max_rbm = 0;
        for (auto i = 0; i < nb_node; ++i)
            if (ranking_scores[i] > max_rbm){
                max_rbm = ranking_scores[i];
            if (ranking_scores[i] < min_cbm)
                min_rbm = ranking_scores[i];
        }
    }

    void MLPH::compute_correlation_based_measure(){

        float mean_y = 0.0;
        for (auto i = 0; i < sample_size; ++i){
            mean_y += objs[i];
        }

        mean_y = mean_y/sample_size;
        std::vector<float> diff_y = std::vector<float>(sample_size);
        float variance_y = 0.0, sum_diff_y = 0.0;
        for (auto i = 0; i < sample_size; ++i){
            diff_y[i] = objs[i] - mean_y;
            variance_y += diff_y[i]*diff_y[i];
            sum_diff_y += diff_y[i];
        }

        std::vector<float> mean_x(nb_node, 0.0);
        std::vector<float> S1(nb_node, 0.0);
        for (auto i = 0; i < sample_size; ++i){
            for (auto j = 0; j < mis_set[i].size(); ++j){
                mean_x[mis_set[i][j]] += 1.0/sample_size;
                S1[mis_set[i][j]] += diff_y[i];
            }
        }

        std::vector<float> variance_x(nb_node, 0.0);
        std::vector<float> variance_xy(nb_node, 0.0);
        corr_xy = std::vector<float>(nb_node, 0.0);
        min_cbm = 1.0; max_cbm = -1.0;

        for (auto i = 0; i < nb_node; ++i){
            variance_x[i] = mean_x[i]*(1-mean_x[i])*sample_size;
            variance_xy[i] = (1-mean_x[i])*S1[i] - mean_x[i]*(sum_diff_y - S1[i]);
                if (mean_x[i] == 0){
                    corr_xy[i] = -1;
                }
                if (mean_x[i] == 1){
                    corr_xy[i] = 1;
                }
                if (variance_x[i] > 0){
                    auto tmp = variance_x[i]*variance_y;
                    if (tmp > 0)
                        corr_xy[i] = variance_xy[i]/sqrt(tmp);
                    else
                        corr_xy[i] = 0;
                }
                if (corr_xy[i] < min_cbm)
                    min_cbm = corr_xy[i];
                if (corr_xy[i] > max_cbm)
                    max_cbm = corr_xy[i];

                if (corr_xy[i]!=corr_xy[i]){
                    std::cout << variance_x[i] << " " << variance_y << "\n";
                }
        }
    }

    void MLPH::make_prediction(){
        compute_correlation_based_measure();
        compute_ranking_based_measure();
        predicted_value = std::vector<float>(nb_node, 0);
        float projection;

        if (max_cbm==0) max_cbm=1;
        
        for (int i = 0; i < nb_node; ++i){

            projection = 
            c0 * ranking_scores[i] / max_rbm +
            c1 * corr_xy[i] / max_cbm + 
            c2 * dual_values[i] / max_dual + 
            c3 * bound_norm[i] + 
            c4 * degree_norm[i] + b;

            predicted_value[i] = 1-1.0/(1+exp(b0*projection+b1));
        }
    }

    void MLPH::run(){
        start_time=get_wall_time();
        random_sampling();
        this->make_prediction();
        this->_run();
        compute_statistics();
    }


    void MLPH::_run(){

        long time_seed = current_time_for_seeding();        
        #pragma omp parallel
        {
            int threadnum = omp_get_thread_num();
            mt19937 mt(time_seed + threadnum*1e5);
            uniform_real_distribution<> dist(0.,1.);
            std::vector<double> distribution (nb_node, 0);

            #pragma omp for
            for(int idx = 0; idx < sample_size; ++idx) {
                if (get_wall_time() - start_time > cutoff){
                    idx = sample_size; continue;
                }

                int nb_candidates = nb_node;
                std::vector<int> candidates(nb_candidates);
                std::iota(candidates.begin(), candidates.end(), 0);
                int j, sel_idx;
                std::vector<int> sample;
                double sum = 0.;
                double obj = 0;;
                double r;
                double prob;
                while (nb_candidates > 0){
                    r = dist(mt);

                    if (nb_candidates == nb_node){
                        sel_idx = idx % nb_node;
                    } else{
                        // calculate sampling distribution;
                        sum = 0.;
                        for (int i = 0; i < nb_candidates; i++){
                            sum += predicted_value[candidates[i]];
                        }
                        if (sum < 1e-8)
                        {
                            for (int i = 0; i < nb_node; i++)
                                distribution[candidates[i]] = 1. / nb_candidates;
                        }else{
                            for (int i = 0; i < nb_node; i++)
                                distribution[candidates[i]] = predicted_value[candidates[i]] / sum;
                        }

                        for (j = 0; j < nb_candidates; j++){
                            prob = distribution[candidates[j]];
                            if (r > prob)
                                r -= prob;
                            else
                                break;
                        }
                        sel_idx = (j==nb_candidates) ? 0 : j;
                    }
                    auto sel_node = candidates[sel_idx];
                    sample.push_back(sel_node);
                    obj += dual_values[sel_node];

                    int num = 0;
                    const std::vector<bool>& neighbors = adj_matrix[sel_node];
                    for (int j = 0; j < nb_candidates; ++j){
                        if (neighbors[candidates[j]] == 0 && j != sel_idx){
                            candidates[num] = candidates[j];
                            num++;
                        }
                    }
                    nb_candidates = num;
                }

                std::sort (sample.begin(), sample.end());  
                std::stringstream ss;
                for (auto k = 0; k < sample.size(); k++){
                    ss << sample[k] << " ";
                }
                std::string identity = ss.str();  
                #pragma omp critical
                {
                    if (best_obj < obj)
                        best_obj = obj;
                    bool duplicate = identites.find(identity) != identites.end();
                    if (!duplicate && 1 - obj < THRESHOLD ){
                            neg_rc_cols.push_back(sample);
                            neg_rc_vals.push_back(1-obj);
                            identites.insert(identity);
                    }
                }
            }
        }

        heur_best_reduced_cost = 1-best_obj;
    }
}