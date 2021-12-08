#include "pricer.h"
#include <stdlib.h>     /* abs */
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

        
        if (column_selection == 0)
            add_partial(basic_cols);
        else if (column_selection == 1)
            add_all(basic_cols);
        else if (column_selection == 5){
            replace_existing(basic_cols, lp_basis);
        }else{
            cout << "ERROR: unrecognized column selection method in include_new_cols()!\n";
            assert(false);
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

    // complexity: O(n)
    void Pricer::add_all(vector<vector<int>>& basic_cols){
        cout << "column selection method: add_all\n";
        basic_cols.insert(basic_cols.end(), neg_rc_cols.begin(), neg_rc_cols.end());
    }

    // complexity bounded by sorting: O(n log(n))
    void Pricer::add_partial(vector<vector<int>>& basic_cols){
        assert(sorted_indices.size()==num_neg_rc_col);
        cout << "column selection method: add_partial\n";
        sorted_indices = sort_indexes_inc(neg_rc_vals);
        for (long i = 0; i < neg_rc_cols.size() && i < upper_col_limit; ++i)
            basic_cols.push_back(neg_rc_cols[sorted_indices[i]]);
    }

    void Pricer::replace_existing(vector<vector<int>>& basic_cols, vector<int>& lp_basis){
        cout << "column selection method: replace_existing\n";

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
};


