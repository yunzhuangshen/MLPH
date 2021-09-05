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
        
        if (method_type == METHOD_TYPE::BP_DEF || method_type==METHOD_TYPE::BP_None ||
            method_type == METHOD_TYPE::GUROBI || method_type == METHOD_TYPE::GUROBI_HEUR) return;
        // cout << "complementary graph #nodes/#edges/density: " 
        // << nb_node<< " "<< nb_edge << " " << nb_edge / ((nb_node*(nb_node-1))/2.)<< "\n";

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
                cout << "ERROR k!= Node_Degree in constructer Heur_pricer_misc.cpp\n";
            // assert(k==Node_Degree[i]);
            Node_Bound[i] += Node_Weight[i];
        }
        
    }

    Heur_pricer_misc::~Heur_pricer_misc(){
        if (method_type == METHOD_TYPE::BP_DEF || method_type==METHOD_TYPE::BP_None ||
            method_type == METHOD_TYPE::GUROBI || method_type == METHOD_TYPE::GUROBI_HEUR) return;
        free(Node_Degree);
        free(Node_Bound);
        free(Node_Weight);
        for (auto i = 0; i < nb_node; i++)
            free(AdjacentList[i]);
        free(AdjacentList);
    }

    void Heur_pricer_misc::solve_mwc_tsm(){

        // solve the reduced graph using TSM
        TSM_SOLVER::TSM solver;
        char const *name2 = "ret";
        solver.tsm("", name2, nb_node, nb_edge, cutoff, AdjacentList, Node_Degree, Node_Weight, Node_Bound);
        objs = solver.sol_objs; sols = solver.sols;
    }

    void Heur_pricer_misc::solve_mwc_lscc(){
        LSCC_SOLVER::LSCC solver;
        solver.lscc(nb_node, nb_edge, cutoff, AdjacentList, Node_Degree, Node_Weight);
        objs = solver.sol_objs; sols = solver.sols;
    }

    void Heur_pricer_misc::solve_mwc_fastwclq(){

        long time_seed = current_time_for_seeding();

        FASTWCLQ_SOLVER::FASTWCLQ solver;
        solver.fastwclq(nb_node, time_seed, nb_node, nb_edge, cutoff, AdjacentList, Node_Degree, Node_Weight, Node_Bound);
        objs.insert(objs.end(), solver.sol_objs.begin(), solver.sol_objs.end());
        sols.insert(sols.end(), solver.sols.begin(), solver.sols.end());
    }

    void Heur_pricer_misc::postprocessing(){
        if (method_type == METHOD_TYPE::BP_DEF || method_type==METHOD_TYPE::BP_None ||
            method_type == METHOD_TYPE::GUROBI || method_type == METHOD_TYPE::GUROBI_HEUR) return;
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
    }

    void Heur_pricer_misc::run(){
        if (cutoff <= 0)
            return;
        switch (method_type)
        {
        case METHOD_TYPE::TSM: solve_mwc_tsm(); break;
        case METHOD_TYPE::LSCC: solve_mwc_lscc(); break;
        case METHOD_TYPE::FASTWCLQ: solve_mwc_fastwclq(); break;
        case METHOD_TYPE::GUROBI: 
        case METHOD_TYPE::GUROBI_HEUR: solve_mwis_gurobi(); break;
        case METHOD_TYPE::BP_DEF: solve_mwis_greedy(); break;
        case METHOD_TYPE::BP_None: break;
        default:
            cout << "not recognizied option for heuristic solver:" << method_type << "\n";
            exit(-1);
        }
        postprocessing();
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


    void Heur_pricer_misc::solve_mwis_gurobi() {
        // try{

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

        model.set(GRB_IntParam_Threads, 1);
        model.set(GRB_IntParam_Presolve, 0);
        if (method_type == METHOD_TYPE::GUROBI_HEUR){
            model.set(GRB_DoubleParam_Heuristics, 0.95);
            model.set(GRB_IntParam_PoolSearchMode, 2);
        }

        // the objective
        GRBLinExpr tot = 0;
        for(int i = 0; i < nb_node; ++i){
            tot += x[i] * dual_values[i];
        }
        model.setObjective(tot,GRB_MAXIMIZE);
        model.update();
        model.optimize();

        auto nsol = model.get(GRB_IntAttr_SolCount);
        // cout << "nsol: " << nsol << "\n";
        for (int i = 0; i < nsol; i++){
            vector<int> sol;
            model.set(GRB_IntParam_SolutionNumber, i);
            auto rc = 1. - model.get(GRB_DoubleAttr_PoolObjVal);
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
        //             }catch(GRBException& e){
        //     cout << "message: " << e.getMessage() << "\n";
        // }
    }
}
