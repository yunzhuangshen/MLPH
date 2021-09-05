#include "pricer.h"
#include <stdlib.h>     /* abs */
#include "fastcluster.h"
#include <list>
#include <random>
#include <bits/stdc++.h>
#include <limits>
#include <cassert>
#include <cmath>
#define pow2(n) ( 1 << (n) )

namespace GCP{
using namespace std;

    long Pricer::current_time_for_seeding(){
        using namespace chrono;
        long now = duration_cast< nanoseconds >(
        system_clock::now().time_since_epoch()).count();
        return now;
    }
    
    vector<long> Pricer::sort_indexes_inc(const vector<double> &v) {
        // initialize original index locations
        vector<long> idx(v.size());
        iota(idx.begin(), idx.end(), 0);
        stable_sort(idx.begin(), idx.end(),
            [&v](long i1, long i2) {return v[i1] < v[i2];});

        return idx;
    }

    vector<long> Pricer::sort_indexes_dec(const vector<double> &v) {
        // initialize original index locations
        vector<long> idx(v.size());
        iota(idx.begin(), idx.end(), 0);
        stable_sort(idx.begin(), idx.end(),
            [&v](long i1, long i2) {return v[i1] > v[i2];});
        return idx;
    }

    double Pricer::get_wall_time(){
        struct timeval time;
        if (gettimeofday(&time,NULL)){
            return 0;
        }
        return (double)time.tv_sec + (double)time.tv_usec * .000001;
    }

    void Pricer::compute_statistics(){
                // calculate statistics
        num_neg_rc_col = neg_rc_vals.size();
        if (num_neg_rc_col == 0){
            best_rc = 0.; mean_rc=0.; median_rc=0.; stdev_rc=0.;
        }else{
            sorted_indices = sort_indexes_inc(neg_rc_vals);
            best_rc = neg_rc_vals[sorted_indices[0]]; 
            median_rc = neg_rc_vals[sorted_indices[long(num_neg_rc_col/2)]];
            mean_rc = 0.;
            for (auto i = 0; i < num_neg_rc_col; ++i)
                mean_rc+=neg_rc_vals[i];
            mean_rc /= num_neg_rc_col;
            for (auto i = 0; i < num_neg_rc_col; ++i)
                stdev_rc +=  (neg_rc_vals[i]-mean_rc) * (neg_rc_vals[i]-mean_rc);
            stdev_rc = sqrt(stdev_rc/num_neg_rc_col);
        }
    }

    void Pricer::include_new_cols(vector<vector<int>>& basic_cols, vector<int>& lp_basis){
        if (method >= 7) assert(column_selection==0);
        double t0 = get_wall_time();

        if (column_selection == 4){
            include_new_cols_random_sampling_all(basic_cols, lp_basis);
        } else if (column_selection == 5){
            include_new_cols_nrc_greedy_all(basic_cols, lp_basis);
        }
        else if (num_neg_rc_col <= upper_col_limit){
            cout << "num neg cols is less than upper limit\n";
            include_new_cols_all(basic_cols);
        }else{
        
            if (column_selection == 0)
                include_new_cols_nrc_greedy(basic_cols, lp_basis);
            else if (column_selection == 1)
                include_new_cols_all(basic_cols);
            else if (column_selection == 2)
                include_new_cols_nrc_sampling(basic_cols);
            else if (column_selection == 3)
                include_new_cols_random_sampling(basic_cols);
            else{
                cout << "ERROR: unrecognized column selection method in include_new_cols()!\n";
                assert(false);
            }
        }

        auto duration = get_wall_time() - t0;
        column_selection_time = duration;
        cout << "column selection time: " << duration << "\n";
    }

    vector<double> Pricer::compute_nrc(vector<vector<int>>& basic_cols){
        vector<double> nrcs(basic_cols.size(), 1);
        for (auto i = 0; i < basic_cols.size(); i++){
            for (auto v : basic_cols[i]){
                nrcs[i] -= dual_values[v];
            }
        }
        return nrcs;
    }

    void Pricer::include_new_cols_nrc_greedy_all(vector<vector<int>>& basic_cols, vector<int>& lp_basis){
        
        cout << "column selection method: include_new_cols_nrc_greedy_all\n";

        int nbasics = 0;
        vector<double> tot_nrcs;
        vector<vector<int>>new_basic_cols;
        vector<vector<int>>tot_cols;

        for (auto i = 0; i < lp_basis.size(); i++){
            if (lp_basis[i]==0) {
                nbasics++;
                new_basic_cols.push_back(basic_cols[i]);
            }else {
                double tmp = 1;
                for (auto v : basic_cols[i])
                    tmp -= dual_values[v];
                tot_nrcs.push_back(tmp);
                tot_cols.push_back(basic_cols[i]);
            }
        }

        tot_cols.insert(tot_cols.end(), neg_rc_cols.begin(), neg_rc_cols.end());
        tot_nrcs.insert(tot_nrcs.end(), neg_rc_vals.begin(), neg_rc_vals.end());

        sorted_indices = sort_indexes_inc(tot_nrcs);

        vector<vector<long>> selected_col_indices(nb_node,vector<long>());

        bool stop = false;
        long tot_added_cols=0;
        for(auto i = 0; i < sorted_indices.size() && !stop; i++){
            auto& col = tot_cols[sorted_indices[i]];
            for (auto node : col){
                if (selected_col_indices[node].size()<10){
                    selected_col_indices[node].push_back(sorted_indices[i]);
                    tot_added_cols++;
                    if (tot_added_cols >= 10 * nb_node) stop = true;
                    break;
                }
            }
        }

        for (auto& indices:selected_col_indices){
            for (auto col_idx : indices){
                new_basic_cols.push_back(tot_cols[col_idx]);
            }
        }

        basic_cols.clear(); 
        basic_cols.insert(basic_cols.end(), new_basic_cols.begin(), new_basic_cols.end());
    }

    // complexity bounded by sorting: O(n log(n))
    void Pricer::include_new_cols_nrc_greedy(vector<vector<int>>& basic_cols, vector<int>& lp_basis){
        assert(sorted_indices.size()==num_neg_rc_col);
        cout << "column selection method: include_new_cols_nrc_greedy\n";

        // replace in place
        // vector<vector<int>> new_basic_cols;
        // for (auto i = 0; i < lp_basis.size(); i++){
        //     if (lp_basis[i]==0) new_basic_cols.push_back(basic_cols[i]);
        // }
        
        // auto pre_size = new_basic_cols.size();
        // for (auto i = 0; i < basic_cols.size() - pre_size && i < neg_rc_cols.size(); ++i){
        //         new_basic_cols.push_back(neg_rc_cols[sorted_indices[i]]);
        // }
        // pre_size = new_basic_cols.size();
        // for (auto i = 0; i < basic_cols.size() -pre_size ; ++i){
        //     new_basic_cols.push_back(basic_cols[i]);
        // }
        // basic_cols.clear();
        // basic_cols.insert(basic_cols.end(), new_basic_cols.begin(), new_basic_cols.end());
        sorted_indices = sort_indexes_inc(neg_rc_vals);
        for (long i = 0; i < neg_rc_cols.size() && i < upper_col_limit; ++i)
            basic_cols.push_back(neg_rc_cols[sorted_indices[i]]);
    }

    // complexity: O(n)
    void Pricer::include_new_cols_all(vector<vector<int>>& basic_cols){
        cout << "column selection method: include_new_cols_all\n";
        basic_cols.insert(basic_cols.end(), neg_rc_cols.begin(), neg_rc_cols.end());
    }

    void Pricer::include_new_cols_random_sampling(vector<vector<int>>& basic_cols){
        cout << "column selection method: include_new_cols_random_sampling\n";
        vector<long> tmp(neg_rc_cols.size()); 
        iota(tmp.begin(), tmp.end(), 0);
        mt19937 mt(current_time_for_seeding());
        std::shuffle(tmp.begin(), tmp.end(), mt);

        for (auto i = 0; i < nb_node; i++)
            basic_cols.push_back(neg_rc_cols[tmp[i]]);
    }

    void Pricer::include_new_cols_random_sampling_all(vector<vector<int>>& basic_cols, vector<int>& lp_basis){
        cout << "column selection method: include_new_cols_random_sampling_all\n";
        int nbasics = 0;
        vector<double> tot_nrcs;
        vector<vector<int>>new_basic_cols;
        vector<vector<int>>tot_cols;

        for (auto i = 0; i < lp_basis.size(); i++){
            if (lp_basis[i]==0) {
                nbasics++;
                new_basic_cols.push_back(basic_cols[i]);
            }else {
                double tmp = 1;
                for (auto v : basic_cols[i])
                    tmp -= dual_values[v];
                tot_nrcs.push_back(tmp);
                tot_cols.push_back(basic_cols[i]);
            }
        }

        tot_cols.insert(tot_cols.end(), neg_rc_cols.begin(), neg_rc_cols.end());
        tot_nrcs.insert(tot_nrcs.end(), neg_rc_vals.begin(), neg_rc_vals.end());

        vector<long> tmp(tot_cols.size()); 
        iota(tmp.begin(), tmp.end(), 0);
        mt19937 mt(current_time_for_seeding());
        std::shuffle(tmp.begin(), tmp.end(), mt);

        auto i = 0;
        while (new_basic_cols.size() < basic_cols.size())
            new_basic_cols.push_back(tot_cols[tmp[i]]);
        
        basic_cols.clear();
        basic_cols.insert(basic_cols.end(), new_basic_cols.begin(), new_basic_cols.end());
    }

    // complexity using binary search: O(n + k * log2(n))
    // https://stackoverflow.com/questions/57599509/c-random-non-repeated-integers-with-weights
    void Pricer::include_new_cols_nrc_sampling(vector<vector<int>>& basic_cols){
        cout << "column selection method: include_new_cols_nrc_sampling\n";
        long time_seed = current_time_for_seeding();        
        mt19937 rng(time_seed);
        
        const vector<double>& weights = neg_rc_vals;
        long rnd_max = weights.size();

        /* determine smallest power of two that is larger than N */
        long tree_levels = ceil(log2((double) rnd_max));

        /* initialize vector with place-holders for perfectly-balanced tree */
        std::vector<double> tree_weights(pow2(tree_levels + 1));

        /* compute sums for the tree leaves at each node */
        long offset = pow2(tree_levels) - 1;
        for (long ix = 0; ix < rnd_max; ix++) {
            tree_weights[ix + offset] = max(0., -weights[ix]);
        }
        for (long ix = pow2(tree_levels+1) - 1; ix > 0; ix--) {
            tree_weights[(ix - 1) / 2] += tree_weights[ix];
        }

        /* sample according to uniform distribution */
        double rnd_subrange, w_left;
        double curr_subrange;
        long curr_ix;
        std::vector<long> sampled_col_indices(upper_col_limit);
        for (long el = 0; el < upper_col_limit; el++) {

            /* go down the tree by drawing a random number and
            checking if it falls in the left or right sub-ranges */
            curr_ix = 0;
            curr_subrange = tree_weights[0];
            for (long lev = 0; lev < tree_levels; lev++) {
                rnd_subrange = std::uniform_real_distribution<double>(0, curr_subrange)(rng);
                w_left = tree_weights[2 * curr_ix + 1];
                curr_ix = 2 * curr_ix + 1 + (rnd_subrange >= w_left);
                curr_subrange = tree_weights[curr_ix];
            }

            /* finally, add element from this iteration */
            sampled_col_indices[el] = curr_ix - offset;

            /* now remove the weight of the chosen element */
            tree_weights[curr_ix] = 0;
            for (long lev = 0; lev < tree_levels; lev++) {
                curr_ix = (curr_ix - 1) / 2;
                tree_weights[curr_ix] =   tree_weights[2 * curr_ix + 1]
                                        + tree_weights[2 * curr_ix + 2];
            }
        }

        for (long i = 0; i < upper_col_limit; ++i){
            // cout << sampled_col_indices[i] << "\n";
            basic_cols.push_back(neg_rc_cols[sampled_col_indices[i]]);
        }
    }

    void Pricer::include_new_cols_hierarchical_clustering(vector<vector<int>>& basic_cols){
        assert(sorted_indices.size()==num_neg_rc_col);
        assert(num_neg_rc_col > upper_col_limit);
        cout << "column selection method: include_new_cols_hierarchical_clustering\n";

        double t0 = get_wall_time();
        long val = (num_neg_rc_col*(num_neg_rc_col-1))/2;
        double* distmat = new double[val];
        long k,i,j;
        for (i=k=0; i<num_neg_rc_col; i++) {
            for (j=i+1; j<num_neg_rc_col; j++) {
                // compute distance between two columns  
                distmat[k] = calc_dist(neg_rc_cols[i], neg_rc_cols[j]);
                k++;
            }
        }
        cout << "time for building distance matrix: " << get_wall_time() - t0 << "\n";

        int* merge = new int[2*(num_neg_rc_col-1)];
        double* height = new double[num_neg_rc_col-1];
        hclust_fast(num_neg_rc_col, distmat, HCLUST_METHOD_SINGLE, merge, height);

        int* labels = new int[num_neg_rc_col];
        // partitioning into nclust clusters
        cutree_k(num_neg_rc_col, merge, upper_col_limit, labels);
        // stop clustering at step with custer distance >= cdist
        // cutree_cdist(num_neg_rc_col, merge, height, cdist, labels);

        std::vector<int> best_col_indices_each_cluster(num_neg_rc_col);
        std::vector<double> best_rc_values_each_cluster(num_neg_rc_col, 1e2);

        for (auto i = 0; i < num_neg_rc_col; i++){
            if(neg_rc_vals[i] < best_rc_values_each_cluster[labels[i]]){
                best_rc_values_each_cluster[labels[i]] = neg_rc_vals[i];
                best_col_indices_each_cluster[labels[i]] = i; 
            }
        }
        
        for (auto i = 0; i < upper_col_limit; i++)
            basic_cols.push_back(neg_rc_cols[best_col_indices_each_cluster[i]]);

        delete[] distmat;
        delete[] merge;
        delete[] height;
        delete[] labels;
    }


    void Pricer::include_new_cols_hierarchical_clustering_all(vector<vector<int>>& basic_cols, vector<int>& lp_basis){
        assert(sorted_indices.size()==num_neg_rc_col);
        assert(num_neg_rc_col > upper_col_limit);
        double t0 = get_wall_time();


        vector<vector<int>> tot_cols;
        long ncluster=basic_cols.size();
        vector<double> tot_nrcs;
        for (auto i = 0; i < lp_basis.size(); i++){
            if (lp_basis[i]==0) ncluster--;
            else {
                tot_cols.push_back(basic_cols[i]);
                double tmp = 1;
                for (auto v : basic_cols[i])
                    tmp -= dual_values[v];
                tot_nrcs.push_back(tmp);
            }
        }
        long ntot_cols = ncluster + num_neg_rc_col;
        
        tot_cols.insert(tot_cols.end(),neg_rc_cols.begin(), neg_rc_cols.end());
        tot_nrcs.insert(tot_nrcs.end(), neg_rc_vals.begin(), neg_rc_vals.end());

        long val = (ntot_cols*(ntot_cols-1))/2;
        double* distmat = new double[val];
        long k,i,j;
        double new_col_dist_max=0;
        double existing_col_dist_min=1e10;
        double new_col_to_existing_col_dist_min=1e10;
        for (i=k=0; i<ntot_cols; i++) {
            for (j=i+1; j<ntot_cols; j++) {
                // compute distance between two columns  
                distmat[k] = calc_dist(tot_cols[i], tot_cols[j]);
                // if (i < nbasics && j < nbasics){
                //     if (distmat[k]> new_col_dist_max)
                //         new_col_dist_max = distmat[k];    
                // }else if (i >= nbasics && j >= nbasics){
                //     if (distmat[k] < existing_col_dist_min)
                //         existing_col_dist_min = distmat[k];    
                // }else if (i < nbasics && j >= nbasics){
                //     if (distmat[k] < new_col_to_existing_col_dist_min)
                //         new_col_to_existing_col_dist_min = distmat[k];    
                // }
                k++;
            }
        }

        // cout << "new_col_dist_max: " << new_col_dist_max << "\n";
        // cout << "existing_col_dist_min: " << existing_col_dist_min << "\n";
        // cout << "new_col_to_existing_col_dist_min: " << new_col_to_existing_col_dist_min << "\n";

        cout << "time for computing the distance matrix: " << get_wall_time() - t0 << "\n";
        int* merge = new int[2*(ntot_cols-1)];
        double* height = new double[ntot_cols-1];
        // hclust_fast(ntot_cols, distmat, HCLUST_METHOD_SINGLE, merge, height);
        hclust_fast(ntot_cols, distmat, HCLUST_METHOD_COMPLETE, merge, height);

        int* labels = new int[ntot_cols];
        // partitioning into nclust clusters
        cutree_k(ntot_cols, merge, ncluster, labels);
        // stop clustering at step with custer distance >= cdist
        // cutree_cdist(ntot_cols, merge, height, cdist, labels);

        //TODO select the columns with the maximum reduced costs and add to basis

        std::vector<int> best_col_indices_each_cluster(ncluster);
        std::vector<double> best_rc_values_each_cluster(ncluster, 1e2);

        // int tmp = -1;
        for (auto i = 0; i < ntot_cols; i++){
            if(tot_nrcs[i] < best_rc_values_each_cluster[labels[i]]){
                best_rc_values_each_cluster[labels[i]] = tot_nrcs[i];
                best_col_indices_each_cluster[labels[i]] = i; 
            }
        }

        vector<vector<int>> new_basic_cols;
        int nselected_basic=0;
        int nselected_new=0;

        for (auto i = 0; i < lp_basis.size(); i++){
            if (lp_basis[i]==0) new_basic_cols.push_back(basic_cols[i]);
        }
        for (auto i = 0; i < ncluster; i++){
            new_basic_cols.push_back(tot_cols[best_col_indices_each_cluster[i]]);

            if (best_col_indices_each_cluster[i] >= ncluster) nselected_new ++;
            else nselected_basic++;
        }

        cout << "nselected_basic/nselected_new: " << nselected_basic << " " << nselected_new << "\n";
        basic_cols.clear();
        basic_cols.insert(basic_cols.end(), new_basic_cols.begin(), new_basic_cols.end());
        auto duration = get_wall_time() - t0;
        column_selection_time = duration;
        delete[] distmat;
        delete[] merge;
        delete[] height;
        delete[] labels;
    }


    inline double Pricer::calc_dist(vector<int>& col1, vector<int>& col2){
        double diff = 0;

        // vector<int> col1_binary(nb_node, 0); vector<int> col2_binary(nb_node, 0);
        // for (auto i = 0; i < col1.size(); i++)
        //     col1_binary[col1[i]] = 1;
        // for (auto i = 0; i < col2.size(); i++)
        //     col2_binary[col2[i]] = 1;            
        // for (auto i = 0; i < nb_node; i++)
        //     diff += abs(col1_binary[i] - col2_binary[i]);

        int k1 = 0, k2 = 0, common = 0;
        while (k1 < col1.size() && k2 < col2.size()){
            if(col1[k1] > col2[k2]){
                k2++;
            } else if (col1[k1] < col2[k2]){
                k1++;
            } else {
                k1++, k2++, common++;
            }
        }
        diff = col1.size() + col2.size() - 2 * common;
        return diff;
    }


    inline double Pricer::calc_dist_kmeans(vector<double>& centroid, vector<int>& point){
        double diff = 0;
        assert(point.size()==nb_node);
        assert(centroid.size()==point.size());
        for (auto i = 0; i < point.size(); i++)
            diff += pow(centroid[i] - point[i], 2);
        return sqrt(diff);
    }

    // void Pricer::include_new_cols_kmeans(vector<vector<int>>& basic_cols){

    //     int k = upper_col_limit, npoints = num_neg_rc_col;
    //     int dimension = nb_node; int niterations = 1e8;
    //     double t0 = get_wall_time();
    //     // binary string representating and dist matrix
    //     vector<vector<int>> points(npoints, vector<int>(nb_node, 0));
    //     for (int i = 0; i < npoints; ++i){
    //         for (auto v : neg_rc_cols[i])
    //             points[i][v] = 1;
    //     }

    //     vector<int> cluster_belong_to(npoints, -1);
    //     vector<double> min_dist_to_centroid(npoints, 1e8);
    //     vector<vector<int>>points_in_cluster(k, vector<int>());

    //     // initialize clusters
    //     vector<vector<double>> centroids(k); 
    //     vector<double> tmp(npoints); 
    //     iota(tmp.begin(), tmp.end(), 0);
    //     mt19937 mt(current_time_for_seeding());
    //     std::shuffle(tmp.begin(), tmp.end(), mt);

    //     tmp.resize(k);
    //     for (auto i=0; i<k;i++){
    //         vector<double> vec(points[tmp[i]].begin(), points[tmp[i]].end());
    //         centroids[i] = vec;
    //     }
        

    //     int iteration = 1; 
    //     for (;iteration < niterations; iteration++){

    //         bool done = true;
            
    //         // Add all points to their nearest cluster
    //         vector<int> cluster_belong_to_pre(cluster_belong_to.begin(),cluster_belong_to.end());
    //         for (int cluster_idx = 0; cluster_idx < k; cluster_idx++){
    //             auto centroid = centroids[cluster_idx];
    //             for (auto i = 0; i < npoints; i++){
    //                 auto dist = calc_dist_kmeans(centroid, points[i]);
    //                 if (dist < min_dist_to_centroid[i]){
    //                     min_dist_to_centroid[i] = dist;
    //                     cluster_belong_to[i] = cluster_idx;
    //                 }
    //             }
    //         }

    //         auto nchanged = 0;
    //         for (auto i = 0; i < npoints; i++){
    //             if (cluster_belong_to[i] != cluster_belong_to_pre[i]){
    //                 nchanged++;
    //             }
    //         }
    //         cout << "current iteration: " << iteration << ", nb changes: " << nchanged << "\n";

    //         // check whether converged
    //         for (auto i = 0; i < npoints; i++){
    //             if (cluster_belong_to[i] != cluster_belong_to_pre[i]){
    //                 done=false;break;
    //             }
    //         }
    //         if (done) break;
            
    //         for (auto i = 0; i < k; i++){
    //             points_in_cluster[i].clear();
    //         }
    //         for (auto i = 0; i < npoints; i++)
    //             points_in_cluster[cluster_belong_to[i]].push_back(i);
            

    //         // Recalculating the center of each cluster
    //         for (auto i = 0; i < k; i++){
    //             fill(centroids[i].begin(), centroids[i].end(), 0.);
    //             for (auto point_idx : points_in_cluster[i]){
    //                 for (auto j = 0; j < dimension; j++)
    //                     centroids[i][j] += points[point_idx][j];
    //             }

    //             for (auto j = 0; j < dimension; j++)
    //                 centroids[i][j] /= points_in_cluster[i].size();
    //         }
    //     }

    //     for (auto& cluster : points_in_cluster){
    //         double min_rc=INF; int best_col_idx = -1;
    //         for (auto col_idx : cluster){
    //             if (neg_rc_vals[col_idx] < min_rc){
    //                 min_rc = neg_rc_vals[col_idx];
    //                 best_col_idx = col_idx;
    //             }
    //         }
    //         basic_cols.push_back(neg_rc_cols[best_col_idx]);
    //     }

    //     auto duration = get_wall_time() - t0;
    //     column_selection_time = duration;
    //     cout << "kmeans niterations: " << iteration << "\n";
    //     cout << "clustering time: " << duration << "\n";

    // }


    // void Pricer::include_new_cols_my_hierarchical_clustering(vector<vector<int>>& basic_cols){
    //     assert(sorted_indices.size()==num_neg_rc_col);
    //     assert(num_neg_rc_col > upper_col_limit);

    //     double t0 = get_wall_time();
    //     vector<vector<double>> distmat(num_neg_rc_col, vector<double>(num_neg_rc_col, 0));
    //     vector<int> dist_to_closest(num_neg_rc_col);
    //     for (auto i=0; i<num_neg_rc_col; i++){ 
    //         auto closest_val = INF;
    //         for (auto j=0; j<num_neg_rc_col; j++){
    //             if (i==j) {
    //                 distmat[i][j] = INF;
    //             }else{
    //                 distmat[i][j] = calc_dist(neg_rc_cols[i], neg_rc_cols[j]);
    //                 if (distmat[i][j] < closest_val){
    //                     closest_val = distmat[i][j];
    //                     dist_to_closest[i] = j;
    //                 }
    //             }
    //         }
    //     }

    //     int cur_ncluster=num_neg_rc_col;
        
    //     vector<bool> candidates(num_neg_rc_col, true);

    //     while(cur_ncluster > upper_col_limit){
            
    //         //find the most similar pair in candidates
    //         double min_dist = INF;
    //         int i1;
    //         for(auto i = 0; i < num_neg_rc_col; i++){
    //             if (candidates[i] && min_dist > distmat[i][dist_to_closest[i]]){
    //                 min_dist = distmat[i][dist_to_closest[i]];
    //                 i1=i;   
    //             }
    //         }
    //         int i2 = dist_to_closest[i1];

    //         // i1 col has better nrc and thus it absorbs i2 col
    //         int removed_col;
    //         if (neg_rc_vals[i1] < neg_rc_vals[i2]) {
    //             candidates[i2] = false;
    //             removed_col = i2;
    //             // update distance matrix in i2th row.
    //             for(auto i = 0; i < num_neg_rc_col; i++){
    //                 distmat[i2][i] = distmat[i1][i];
    //                 distmat[i][i2] = distmat[i1][i];
    //             }
    //             distmat[i2][i2] = INFINITY;
    //         }
    //         else{
    //              candidates[i1] = false;
    //             removed_col = i1;
    //             // update distance matrix in i1th row.
    //             for(auto i = 0; i < num_neg_rc_col; i++){
    //                 distmat[i1][i] = distmat[i2][i];
    //                 distmat[i][i1] = distmat[i2][i];
    //             }
    //             distmat[i1][i1] = INFINITY;
    //         }
            
    //         //update d whose closed one to removed_col
    //         for (auto i = 0; i < num_neg_rc_col; i++){
    //             if (candidates[i] && dist_to_closest[i] == removed_col){
    //                 auto closest_val = INF;
    //                 for (auto j=0; j<num_neg_rc_col; j++){
    //                     if (candidates[j] && distmat[i][j] < closest_val){
    //                         closest_val = distmat[i][j];
    //                         dist_to_closest[i] = j;
    //                     } 
    //                 }
    //             }
    //         }

    //         --cur_ncluster;
    //     }

    //     int ctr=0;
    //     for (auto i = 0; i != num_neg_rc_col; i++){
    //         if (candidates[i]){
    //             basic_cols.push_back(neg_rc_cols[i]);
    //             ctr++;
    //         }
    //     }
    //     assert(ctr == upper_col_limit);
    //     auto duration = get_wall_time() - t0;
    //     column_selection_time = duration;
    //     cout << "clustering time: " << column_selection_time << "\n";
    // }

};


