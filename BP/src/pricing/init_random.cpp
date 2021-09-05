#include "init_random.h"
#include <random>

namespace GCP {
    using namespace std;

    Random_col::Random_col(MWISP_INST& inst, int seed, int nb_mwis):
    mis_adj_matrix{inst.adj_matrix}, nb_node{inst.nb_node}, nb_mwis{nb_mwis}
    {
    }

    void Random_col::add_random_cols(SCIP* scip, SCIP_CONS** constraints){

        mt19937 mt(seed);
        uniform_int_distribution<int> dist(0,RAND_MAX);
        SCIP_VAR*        var;                   /* pointer to the new created variable */
        int v, idx, num, setnumber;
        vector<int> candidates(nb_node);
        int nb_candidates;
        vector<int> col;
        for(int i = 0; i < nb_mwis; ++i) {
            col.clear();
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
                col.push_back(v);
                num = 0;
                for (int j = 0; j < nb_candidates; ++j){
                    if (mis_adj_matrix[v][candidates[j]] == 0 && j != idx){
                        candidates[num] = candidates[j];
                        num++;
                    }
                }
                nb_candidates = num;
            }
            sort(col.begin(), col.end(), greater<int>());

            COLORprobAddNewStableSet(scip, col.data(), col.size(), &setnumber);
            assert(setnumber != -1);

            /* create variable for the stable set and add it to SCIP */
            SCIPcreateVar(scip, &var, NULL, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY,
                TRUE, FALSE, NULL, NULL, NULL, NULL, (SCIP_VARDATA*)(size_t)setnumber); /*lint !e571*/

            COLORprobAddVarForStableSet(scip, setnumber, var);
            SCIPaddVar(scip, var);
            SCIPchgVarUbLazy(scip, var, 1.0);

            for(int j = 0; j < col.size(); j++ )
            {
                /* add variable to node constraints of nodes in the set */
                SCIPaddCoefSetppc(scip, constraints[col[j]], var);
            }
        }
    }

}

