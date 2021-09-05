#include "ACO.h"
#include <cmath>
#include <limits>
#include <sys/time.h>
#include <cassert>
#include <math.h> 
#include <sstream>

namespace GCP {

    ACO::ACO(double _cutoff, MWISP_INST& inst): 
        degree_norm{inst.degree_norm}, adj_matrix{inst.adj_matrix}{
        dual_values = inst.dual_values;
        nb_node = inst.nb_node;
        nb_edge = inst.nb_edge;
        cutoff = _cutoff;
        sample_size = 30;
        niterations = floor(nb_node*5/3);
        // sample_size=nb_node;
        compute_bound();
        max_dual = 0.;
        for (auto i = 0; i < nb_node; i++){
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
    }

    void ACO::run_iteration(int ith_iteration){
        long time_seed = current_time_for_seeding();    
        vector<double> delta_tau(nb_node, 0.0);
    
        mt19937 mt(time_seed);
        uniform_real_distribution<> dist(0.,1.);
        std::vector<double> distribution (nb_node, 0);
        
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
                    identites.insert(identity);
                    neg_rc_cols.push_back(sample);
                    neg_rc_vals.push_back(1-obj);
                }
            }
            if (best_obj < obj)
                best_obj = obj;

            // update delta_tau_local
            for (auto k = 0u; k < sample.size(); ++k){
                delta_tau[sample[k]] += 1/(1+1000 * (best_obj-obj));
            }
            
        }
        
        // update tau
        for (auto k = 0u; k < nb_node; ++k){
            tau[k] = rho * tau[k] + delta_tau[k];
        }
    }
}
