#include "pricer.h"
#include <stdlib.h>     /* abs */
#include <list>
#include <random>
#include <bits/stdc++.h>
#include <limits>
#include <cassert>
#include <cmath>
#include <bits/stdc++.h>

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

    void Pricer::add_columns(SCIP_PRICERDATA* pricerdata, bool is_root){

        SCIP_VAR*        var;                   /* pointer to the new created variable */
        int              setnumber;             /* index of the new created variable */
        int nb_node = COLORprobGetNNodes(pricerdata->scip);
        double factor = is_root ? pricerdata->root_column_limit_factor : pricerdata->child_column_limit_factor;
        int nb_cols = min((int)neg_rc_cols.size(), int(factor * nb_node));
        std::vector<long> sorted_indices = sort_indexes_inc(neg_rc_vals);
        pricerdata->improving_mwiss.resize(nb_cols);
        best_rc = neg_rc_vals[sorted_indices[0]];
        for (int i = 0; i < nb_cols; i++ )
        {
            pricerdata->improving_mwiss[i] = neg_rc_cols[sorted_indices[i]];
            sort(pricerdata->improving_mwiss[i].begin(), pricerdata->improving_mwiss[i].end(), greater<int>());

            int cur_mwis_nnodes = pricerdata->improving_mwiss[i].size();
            if ( cur_mwis_nnodes > 0 )
            {
                // std::cout << "neg_rc_vals: " << neg_rc_vals[sorted_indices[i]] << "\n";
                if ( SCIPisFeasLT(pricerdata->scip, neg_rc_vals[sorted_indices[i]], 0.0) )
                {
                    int setnumber = -1;
                    /* insert new variable */
                    COLORprobAddNewStableSet(pricerdata->scip, pricerdata->improving_mwiss[i].data(),
                        pricerdata->improving_mwiss[i].size(), &setnumber);
                    /* only insert, if there yet is no variable for this stable set */
                    // std::cout << "setnumber: " << setnumber << "\n";
                    if ( setnumber >= 0  )
                    {
                        /* create variable for the stable set and add it to SCIP */
                        SCIPcreateVar(pricerdata->scip, &var, NULL, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY,
                                TRUE, TRUE, NULL, NULL, NULL, NULL, (SCIP_VARDATA*)(size_t)setnumber);

                        COLORprobAddVarForStableSet(pricerdata->scip, setnumber, var);
                        // SCIPvarMarkDeletable(var);
                        SCIPaddPricedVar(pricerdata->scip, var, 1.0);
                        SCIPchgVarUbLazy(pricerdata->scip, var, 1.0);

                        /* add variable to the constraints in which it appears */
                        for (int j = 0; j < pricerdata->improving_mwiss[i].size(); j++ )
                        {
                            /* add variable to node constraints of nodes in the set */
                            SCIPaddCoefSetppc(pricerdata->scip,
                                    pricerdata->constraints[pricerdata->improving_mwiss[i][j]], var);
                        }
                    }
                }
            }
        }


    }
};


