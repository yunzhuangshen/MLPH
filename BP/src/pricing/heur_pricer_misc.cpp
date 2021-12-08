#include "heur_pricer_misc.h"

using namespace std;

// #define DEBUG 

namespace GCP {
    using namespace std;

    Heur_pricer_misc::Heur_pricer_misc(METHOD_TYPE _method_type, double _cutoff, MWISP_INST& inst):
    mis_adj_list{inst.adj_list}, mis_adj_matrix{inst.adj_matrix}
    {
        dual_values = inst.dual_values;
        method_type = _method_type;
        cutoff=_cutoff;
        nb_node = dual_values.size();
        nb_edge = (nb_node*(nb_node-1))/2. - inst.nb_edge;
    }

    void Heur_pricer_misc::run(){
        if (cutoff <= 0)
            return;
        switch (method_type)
        {
        case METHOD_TYPE::BP_DEF: solve_mwis_greedy(); break;
        case METHOD_TYPE::BP_None: break;
        default:
            cout << "not recognizied option for heuristic solver:" << method_type << "\n";
            exit(-1);
        }
    }

    void Heur_pricer_misc::solve_mwis_greedy(){

        vector<long> sortednodes = sort_indexes_dec(dual_values);
        vector<int> mwis;
        double mwis_obj;
        bool indnode;
        /* insert first node */
        mwis.push_back(sortednodes[0]);
        mwis_obj = dual_values[sortednodes[0]];

        for (int i = 1; i < nb_node; i++ )
        {
            /* test if node is independant to nodes in stable set */
            indnode = true;
            for (int j = 0; j < mwis.size(); j++ )
            {
                if ( mis_adj_matrix[sortednodes[i]][mwis[j]] )
                {
                    indnode = false;
                    break;
                }
            }
            /* if node is independant to nodes in stable set, insert it into stable set*/
            if ( indnode )
            {
                mwis.push_back(sortednodes[i]);
                mwis_obj += dual_values[sortednodes[i]];
            }
        }
        double rc = 1 - mwis_obj;
        if (rc < THRESHOLD){
            neg_rc_cols.push_back(mwis);
            neg_rc_vals.push_back(rc);
        }
    }
}
