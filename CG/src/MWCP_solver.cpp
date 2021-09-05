#include "MWCP_solver.h"
#include <cmath>
#include <limits>
#include <numeric>      // iota
#include <cassert>
#include <set>
#include "gurobi_c++.h"
#include <omp.h>

using namespace std;

// #define DEBUG 

namespace GCP {
    using namespace std;

    MWCP_solver::MWCP_solver(int _method, double _cutoff, const vector<double>& _dual_values, const vector<vector<int>>& _mis_adj_list, const vector<vector<bool>>& _mis_adj_matrix, long long nb_edges, long long _upper_col_limit):
    mis_adj_list{_mis_adj_list}, mis_adj_matrix{_mis_adj_matrix}
    {
        dual_values = _dual_values;
        method = _method;
        cutoff=_cutoff;
        upper_col_limit=_upper_col_limit;
        nb_node = dual_values.size();
        nb_edge = (nb_node*(nb_node-1))/2. - nb_edges;

        if (method==11 || method == 12 || method == 13) return;
        
        cout << "complementary graph #nodes/#edges/density: " 
        << nb_node<< " "<< nb_edge << " " << nb_edge / ((nb_node*(nb_node-1))/2.)<< "\n";

        // convert objective from double to long long
        Node_Weight = (long long *) malloc((nb_node) * sizeof(long long));
        for (auto i = 0; i < nb_node; ++i){
            Node_Weight[i] = dual_values[i] * 1e12;
            // if (Node_Weight[i]==0L) Node_Weight[i]=1L;
        }
        // cout << "size of: " << sizeof(long long) << " "  <<sizeof(double) << "\n";
        // calculate complementary graph 
        Node_Degree = (long long *) malloc((nb_node) * sizeof(long long));
        AdjacentList = (long long **) malloc((nb_node) * sizeof(long long *));
        Node_Bound = (long long *) malloc((nb_node) * sizeof(long long));

        for (auto i = 0; i < nb_node; i++){
            Node_Degree[i] = nb_node - mis_adj_list[i].size() - 1;
            AdjacentList[i] = (long long *) malloc(Node_Degree[i] * sizeof(long long));
            vector<bool> candidates(nb_node, true);
            candidates[i]=false;
            for (auto j = 0; j < mis_adj_list[i].size(); j++){
                candidates[mis_adj_list[i][j]] = false;
            }
            auto k = 0;
            for (auto j = 0; j < nb_node; j++){
                if(candidates[j]){
                    AdjacentList[i][k++] = j;
                    Node_Bound[i] += Node_Weight[j];
                }
            }
            // std::cout << i << " " << k << " " << Node_Degree[i] << "\n";
            if (k!=Node_Degree[i]) 
                cout << "ERROR k!= Node_Degree in constructer MWCP_solver.cpp\n";
            // assert(k==Node_Degree[i]);
            Node_Bound[i] += Node_Weight[i];
        }
        
    }

    MWCP_solver::~MWCP_solver(){
        if (method == 11 || method == 12 || method == 13) return;
        free(Node_Degree);
        free(Node_Bound);
        free(Node_Weight);
        for (auto i = 0; i < nb_node; i++)
            free(AdjacentList[i]);
        free(AdjacentList);
    }

    void MWCP_solver::solve_mwc_tsm(){

        // solve the reduced graph using TSM
        TSM::TSM solver;
        char const *name2 = "ret";
        best_obj = solver.tsm("", name2, nb_node, nb_edge, cutoff, AdjacentList, Node_Degree, Node_Weight, Node_Bound);
        objs = solver.sol_objs; sols = solver.sols;
        isOptimal = solver.isOptimal;
        exact_rc = 1 - best_obj/1e12;
    }

    void MWCP_solver::solve_mwc_lscc(){
        LSCC::LSCC solver;
        best_obj = solver.lscc(nb_node, nb_edge, cutoff, AdjacentList, Node_Degree, Node_Weight);
        objs = solver.sol_objs; sols = solver.sols;
    }

    void MWCP_solver::solve_mwc_fastwclq(){
        long thread_limit = omp_get_max_threads();
        if (thread_limit > 1)
            cout << "fastwclq is paralleled using " << thread_limit << " threads\n";
        long ntries_per_thread = 50. * nb_node / thread_limit;
        long time_seed = current_time_for_seeding();

        #pragma omp parallel for
        for (auto t = 0; t < thread_limit; t++){
            FASTWCLQ::FASTWCLQ solver;
            best_obj = solver.fastwclq(ntries_per_thread, time_seed + t * 1e8, nb_node, nb_edge, cutoff, AdjacentList, Node_Degree, Node_Weight, Node_Bound);
            #pragma omp critical
            {
                objs.insert(objs.end(), solver.sol_objs.begin(), solver.sol_objs.end());
                sols.insert(sols.end(), solver.sols.begin(), solver.sols.end());
            }
        }
    }

    void MWCP_solver::postprocessing(){
        if (method==11 || method == 12 || method == 13) return;

        assert(objs.size()==sols.size());
        auto size_before = objs.size();
        set<string> identities;
        vector<double> objs_distinct; 
        vector<vector<int>> sols_distinct;

        // remove duplicates
        for (auto i = 0; i < size_before; i++){
            stable_sort(sols[i].begin(), sols[i].end());  
            stringstream ss;
            for (auto k = 0; k < sols[i].size(); k++)
                ss << sols[i][k] << " ";
            
            string identity = ss.str();                    
            if (identities.find(identity) == identities.end()){
                objs_distinct.push_back(objs[i]);
                vector<int> vec;
                for (auto v : sols[i])
                    vec.push_back(v);
                sols_distinct.push_back(vec);
            }
            identities.insert(identity);
        }
        
        // vertices in a solution should be independent
        vector<bool> isValid(sols_distinct.size(), true); 
        for (auto k = 0; k < sols_distinct.size(); k++){
            auto& sol = sols_distinct[k];
            for (auto i = 0; i < sol.size();i++){
                auto adj = mis_adj_list[sol[i]];
                for (auto j = i+1;j < sol.size(); j++){
                    if (find(adj.begin(), 
                            adj.end(), sol[j])!=adj.end()){
                        isValid[k] = false;
                        // cout << i+1 << " " << j+1 << "\n";
                        // cout << "ERROR: invalid solution by method " << method << "\n" ;
                    }
                    // assert( find(adj.begin(), 
                    //         adj.end(), sol[j])==adj.end());
                }
            }
        }

        sols.clear(); objs.clear();
        for (auto i = 0; i < sols_distinct.size(); i++){
            if (isValid[i]){
                double rc = 1;
                for (auto v : sols_distinct[i]){
                    rc -= dual_values[v];
                }

                if (rc < THRESHOLD){
                    neg_rc_vals.push_back(rc);
                    neg_rc_cols.push_back(sols_distinct[i]);
                }
            }
        }

        num_neg_rc_col  = neg_rc_vals.size();
        if (method==10 && num_neg_rc_col > 0 && isOptimal){
            int best_col=-1; double best_rc = 1e2;
            for (auto i = 0; i < num_neg_rc_col; i++){
                if (neg_rc_vals[i] < best_rc){
                    best_rc = neg_rc_vals[i];
                    best_col = i;
                }
            }
            optimal_mis = neg_rc_cols[best_col];
        }
    }

    void MWCP_solver::run(){
        if (cutoff <= 0)
            return;
        switch (method)
        {
        case 8: solve_mwc_tsm(); break;
        case 9: solve_mwc_lscc(); break;
        case 10: solve_mwc_fastwclq(); break;
        case 11: 
        case 12: solve_mwis_gurobi(); break;
        case 13: solve_mwis_greedy(); break;
        default:
            cout << "not recognizied option for specialized solver\n";
            exit(-1);
        }
        postprocessing();
        compute_statistics();
    }

    void MWCP_solver::solve_mwis_greedy(){

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

    void MWCP_solver::solve_mwis_gurobi() {

        // setup the model now
        GRBEnv *env;
        vector<GRBVar> x;
        env = new GRBEnv();
        GRBModel model = GRBModel(*env);
        model.set(GRB_IntParam_PoolSolutions, 1e8); //retain all solutions
        model.set(GRB_DoubleParam_TimeLimit, cutoff);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        model.set(GRB_StringAttr_ModelName, "MIP_MWIP");
        // Create variables and set them to be binary
        x.resize(nb_node);
        for (int i = 0; i < nb_node; ++i){
            x[i] = model.addVar(0,1,0,GRB_BINARY);
        }
        model.update();

        // adjacent vertices cannot be selected simultaneously.
        for (int i = 0; i < nb_node; ++i){
            for (int j = i+1; j < nb_node; ++j){
                if (mis_adj_matrix[i][j] == 1){
                    model.addConstr(x[i] + x[j] <= 1, "");
                }
            }
        }
        model.update();
        cout << omp_get_max_threads() << "\n";

        model.set(GRB_IntParam_Threads, omp_get_max_threads());
        model.set(GRB_IntParam_Presolve, 0);
        if (method == 12){
            model.set(GRB_DoubleParam_Heuristics, 0.95);
            model.set(GRB_IntParam_PoolSearchMode, 2);
        }

        // the objective
        GRBLinExpr tot = 1;
        for(int i = 0; i < nb_node; ++i){
            tot -= x[i] * dual_values[i];
        }
        model.setObjective(tot,GRB_MINIMIZE);
        model.update();
        model.optimize();

        heur_best_reduced_cost = model.get(GRB_DoubleAttr_ObjVal);
        isOptimal = (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL);
        exact_rc = model.get(GRB_DoubleAttr_ObjVal);

        auto nsol = model.get(GRB_IntAttr_SolCount);
        // cout << "nsol: " << nsol << "\n";
        for (int i = 0; i < nsol; i++){
            vector<int> sol;
            model.set(GRB_IntParam_SolutionNumber, i);
            auto rc = model.get(GRB_DoubleAttr_PoolObjVal);
            // std:: cout << rc << "\n";
            if (rc < THRESHOLD){
                neg_rc_vals.push_back(rc);
                for (int i = 0; i < nb_node; ++i){
                    if (x[i].get(GRB_DoubleAttr_Xn) > 0.5){
                        sol.push_back(i);
                    }
                }
                neg_rc_cols.push_back(sol);
            }
        }
        
        delete env;
    }
    
    // double MWCP_solver::solve_mwc_wlmc(){
    //     best_obj = wlmc(nb_node, nb_edge, cutoff, AdjacentList, Node_Degree, Node_Weight, Node_Bound);
    //     return best_obj;
    // }

}
