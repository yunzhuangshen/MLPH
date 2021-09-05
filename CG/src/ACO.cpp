#include "ACO.h"
#include <cmath>
#include <limits>
#include <sys/time.h>
#include <cassert>
#include <math.h> 
#include <sstream>
#include <omp.h>

namespace GCP {

    ACO::ACO(int _method, double _cutoff, int _n, int _nedge, int _sample_size, const vector<double>& _degree, 
                const vector<vector<bool>>& _adj_matrix, const vector<double>& _dual_values,
                int _upper_col_limit): 
        degree{_degree}, adj_matrix{_adj_matrix}{
        dual_values = _dual_values;
        method = _method;
        nb_node = _n;
        nb_edge = _nedge;
        cutoff = _cutoff;
        sample_size = 30;
        upper_col_limit = _upper_col_limit;
        niterations = floor(nb_node*5/3);
        cout << "aco iteration: " << niterations << "\n";
        sample_size=nb_node;
        best_rc_current_iteration = vector<double>(niterations);
        num_neg_rc_current_iteration = vector<long>(niterations);
        heur_best_reduced_cost = 1;
        compute_bound();
        max_dual = 0.;
        for (auto i = 0; i < nb_node; i++){
            // cout << i << " " << dual_values[i] << "\n";
            if (dual_values[i] > max_dual)
                max_dual = dual_values[i];
        }
    }

    void ACO::compute_bound(){
        bound = std::vector<double> (nb_node, 0);
        for (auto i = 0; i < nb_node; i++){
            for (auto j = 0; j < nb_node; j++){
                if (i!=j && adj_matrix[i][j] == 0)
                    bound[i] += dual_values[j]; 
            }
        }

        bound_sum = 0;
        for (auto i = 0; i < nb_node; i++)
            bound_sum+=dual_values[i];
    }

    void ACO::run(){
        start_time=get_wall_time();
        objs = vector<double> (sample_size);
        mis_set = vector<vector<int>> (sample_size);
        tau = vector<double>(nb_node, 1.);
        for (auto i = 0; i < nb_node; i++)
            tau[i] += dual_values[i];
        // tau = vector<double>(nb_node, 0.1);
        // for (auto i = 0; i < nb_node; i++)
        //     tau[i] += dual_values[i];

        for (auto i = 0; i < niterations; ++i){
            this->run_iteration(i);
            num_neg_rc_col = 0;
            for (auto idx = 0; idx < sample_size; idx++){
                if (1 - objs[idx] < THRESHOLD){
                    num_neg_rc_col++;
                }
            }
            // cout << "iter: " << i << ", best obj: " << best_obj << "\n";
            best_rc_current_iteration[i] = heur_best_reduced_cost;
            num_neg_rc_current_iteration[i] = num_neg_rc_col;
            if (get_wall_time() - start_time > cutoff)
                break;

            // update dynamic parameters
            if (i < 100) alpha = 1;
            else if (i < 400) alpha = 2;
            else if (i < 800) alpha = 3;
            else alpha = 4;

            T = (1-beta) * T;
            if (rho > 0.95)
                rho = (1-delta_rho) * rho;
            else
                rho = 0.95;
        }
        
        compute_statistics();
    }

    void ACO::run_iteration(int ith_iteration){
        long time_seed = current_time_for_seeding();    
        vector<double> delta_tau(nb_node, 0.0);
    
        #pragma omp parallel
        {
            int threadnum = omp_get_thread_num();
            mt19937 mt(time_seed + threadnum*1e5);
            uniform_real_distribution<> dist(0.,1.);
            std::vector<double> distribution (nb_node, 0);
            
            #pragma omp for
            for (auto idx = 0u; idx < sample_size; ++idx){
                if (get_wall_time() - start_time > cutoff){
                    idx=sample_size;continue;
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

                    if (nb_candidates == nb_node)
                        sel_idx = idx % nb_node;
                    else{
                        // calculate sampling distribution;
                        sum = 0.;
                        for (int i = 0; i < nb_candidates; i++){
                            sum += pow(tau[candidates[i]], alpha);
                        }
                        if (sum < 1e-8)
                        {
                            for (int i = 0; i < nb_node; i++)
                                distribution[candidates[i]] = 1. / nb_candidates;
                        }else{
                            for (int i = 0; i < nb_node; i++)
                                distribution[candidates[i]] = pow(tau[candidates[i]], alpha) / sum;
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
                    // cout << sel_node << " " << dual_values[sel_node] << "\n";
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

                objs[idx] = obj;
                mis_set[idx].resize(sample.size());
                std::copy(sample.begin(), sample.end(), mis_set[idx].begin());

                if (1 - obj < THRESHOLD){
                    std::sort (sample.begin(), sample.end());  
                    std::stringstream ss;
                    for (auto k = 0; k < sample.size(); k++){
                        ss << sample[k] << " ";
                    }
                    std::string identity = ss.str();   

                    if (identites.find(identity) == identites.end()){
                        #pragma omp critical
                        {
                            identites.insert(identity);
                            neg_rc_cols.push_back(sample);
                            neg_rc_vals.push_back(1-obj);
                        }
                    }
                }


                #pragma omp critical
                {
                    if (best_obj < obj)
                        best_obj = obj;

                    // update delta_tau_local
                    for (auto k = 0u; k < sample.size(); ++k){
                        delta_tau[sample[k]] += 1/(1+1000 * (best_obj-obj));
                    }
                }
            }
        }
        
        // update tau
        for (auto k = 0u; k < nb_node; ++k){
            tau[k] = rho * tau[k] + delta_tau[k];
        }
        heur_best_reduced_cost = 1-best_obj;
    }

}
