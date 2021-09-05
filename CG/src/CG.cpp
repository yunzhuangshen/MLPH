#include "CG.h"
#include <cmath>
#include <limits>
#include "MLPH.h"
#include <cassert>
#include <iostream>
#include "pricer.h"
#include "MWCP_solver.h"
#include "ACO.h"
#include <omp.h>
#include <chrono>
#include <cmath>

namespace GCP {

    static double get_wall_time(){
        struct timeval time;
        if (gettimeofday(&time,NULL)){
            return 0;
        }
        return (double)time.tv_sec + (double)time.tv_usec * .000001;
    }

    static long current_time_for_seeding(){
        using namespace chrono;
        long now = duration_cast< nanoseconds >(
        system_clock::now().time_since_epoch()).count();
        return now;
    }

    CG::CG(const Instance& g, double cutoff, double pricer_cutoff, double _thread_limit, int _seed) : g{g}, cutoff{cutoff}, pricer_cutoff{pricer_cutoff}, seed{_seed} {
        thread_limit=_thread_limit;
        omp_set_num_threads(thread_limit);
        initializing_parameters();
    }

    void CG::initializing_parameters(){
        num_mis = mis_factor * g.size();
        cout << "number of MISs generated is " <<num_mis << endl;
        mis_set.resize(num_mis);
        adj_list = g.get_adj_list();
        adj_matrix = vector<vector<bool>>(g.size(), vector<bool>(g.size(), 0));
        for (int i = 0; i < g.size(); ++i){
            for (int j = 0; j < adj_list[i].size(); ++j){
                adj_matrix[i][adj_list[i][j]] = 1;
                adj_matrix[adj_list[i][j]][i] = 1;
            }
        }
        dual_value.resize(g.size());
    }

    void CG::initializing_mis() {

        
        #pragma omp parallel
        {
            int threadnum = omp_get_thread_num();
            mt19937 mt(seed + threadnum*1e5);
            uniform_int_distribution<int> dist(0,RAND_MAX);

            int v, idx, num;
            vector<int> candidates(g.size());
            int nb_candidates;

            #pragma omp for
            for(int i = 0; i < num_mis; ++i) {
                nb_candidates = g.size();
                for (int j = 0; j < nb_candidates; ++j){
                    candidates[j] = j;
                }
                while (nb_candidates > 0){
                    if (nb_candidates == g.size()){
                        idx = i % g.size();
                    } else{
                        idx = dist(mt) % nb_candidates;
                    }
                    v = candidates[idx];
                    mis_set[i].push_back(v);
                    num = 0;
                    for (int j = 0; j < nb_candidates; ++j){
                        if (adj_matrix[v][candidates[j]] == 0 && j != idx){
                            candidates[num] = candidates[j];
                            num++;
                        }
                    }
                    nb_candidates = num;
                }
            }
        }
        
    }

    void CG::solve_restricted_master_problem(bool warm_start) {
//        cout << "computing restricted master problem " << endl;
//        vector<vector<int>> mis_set = g.get_mis_set();
        num_mis = mis_set.size();
        vector<vector<bool>> mis_set_binary(num_mis, vector<bool>(g.size(), 0));
        long v;
        for (long i = 0; i < num_mis; ++i){
            for (long j = 0; j < mis_set[i].size(); ++j){
                v = mis_set[i][j];
                mis_set_binary[i][v] = 1;
            }
        }
        // setup the model now
        try{
            GRBEnv *env;
            vector<GRBVar> x;
            env = new GRBEnv();
            GRBModel model = GRBModel(*env);
            model.getEnv().set(GRB_IntParam_OutputFlag, 0);
            model.set(GRB_IntParam_Threads, thread_limit);
            model.set(GRB_StringAttr_ModelName, "RMP_GCP");
            // Create variables and set them to be binary
            x.resize(num_mis);
            for (long i = 0; i < num_mis; ++i){
                x[i] = model.addVar(0,1,0,GRB_CONTINUOUS);
            }
            model.update();

            // each vertex is covered by at least one set
            vector<GRBConstr> y;
            y.resize(g.size());
            for (long j = 0; j < g.size(); ++j){
                GRBLinExpr rtot = 0;
                for (long i = 0; i < num_mis; ++i){
                    rtot += mis_set_binary[i][j] * x[i];
                }
                y[j] = model.addConstr(rtot >= 1, "");
            }

            if (warm_start)
            {
                auto tmp = lp_vbasis.size();
                for (auto i = 0; i < num_mis - tmp; i++)
                    lp_vbasis.push_back(-1);
                assert(num_mis==lp_vbasis.size());
                for (auto i = 0; i < num_mis; ++i)
                    x[i].set(GRB_IntAttr_VBasis, lp_vbasis[i]);
                tmp = lp_cbasis.size();
                for (auto i = 0; i < g.size() - tmp; i++)
                    lp_cbasis.push_back(-1);
                for (auto j = 0; j < g.size(); ++j)
                    y[j].set(GRB_IntAttr_CBasis, lp_cbasis[j]);
            }
            model.update();
            model.set(GRB_IntParam_Presolve, 0);
            // the objective
            GRBLinExpr tot=0;
            for(long i = 0; i < num_mis; ++i){
                tot += x[i];
            }
            model.setObjective(tot,GRB_MINIMIZE);
            model.update();


            model.optimize();
            lp_bound = model.get(GRB_DoubleAttr_ObjVal);
            cout << "LP bound is " << lp_bound << endl;

            // get optimal dual value
            lp_cbasis.resize(g.size());
            for (long j = 0; j < g.size(); ++j){
                dual_value[j] = y[j].get(GRB_DoubleAttr_Pi);
                lp_cbasis[j] = y[j].get(GRB_IntAttr_CBasis);
            }

            long k = 0;
            lp_vbasis.resize(num_mis);
            for (long i = 0; i < num_mis; ++i){
                lp_vbasis[i] = x[i].get(GRB_IntAttr_VBasis);
                k+=(lp_vbasis[i]==0);
            }
            cout << "lp solution # basic variables: " << k << "/" << lp_vbasis.size() << "\n";
            delete env;
        }catch(GRBException e){
            std::cout << "Gurobi Exception\n";
            cout << e.getErrorCode() << " " << e.getMessage() << "\n";
        }
    }

    bool CG::solve_mwis_tsm(double cutoff, double& min_rc){
        if (cutoff <= 0)
            return false;
        MWCP_solver tsm(8, cutoff, dual_value, adj_list, adj_matrix, g.get_nb_edges(),1e8);
        tsm.run();
        optimal_mis=tsm.optimal_mis;
        min_rc=tsm.exact_rc;
        
        return tsm.isOptimal;
    }

    void CG::collect_training_data(vector<vector<double>>& obj_coef, vector<vector<bool>>& solution){
        initializing_mis();

        min_reduced_cost_exact = -1.0;

        while(min_reduced_cost_exact < -0.000001){
            
            cout << "solving restricted master problem\n";
            solve_restricted_master_problem(false);
            cout << "solving pricing problem\n";
            solve_mwis_tsm(1e8, min_reduced_cost_exact);
            // record training data
            if (cg_iters % 5==0){

                vector<bool> opt_sol(g.size(), false);
                for (auto v : optimal_mis){
                    opt_sol[v]=true;
                }

                obj_coef.push_back(dual_value);
                solution.push_back(opt_sol);
            }

            // add new columns
            mis_set.push_back(optimal_mis);
            cout << "minimum reduced cost is " << min_reduced_cost_exact << endl;

            if (cg_iters++>=25) break;
        }
    }


    double CG::optimize_LM(int method, double b0, double b1){
        double start_time = get_wall_time();

        initializing_mis();
        min_reduced_cost_exact = -1.0;

        double t0;
        double duration;
        t0 = get_wall_time();
        solve_restricted_master_problem(false);
        duration = get_wall_time()-t0;
        time_duration_master += duration;
        cout << "time used: " << duration << "\n";
        cout << "solve pricing mwis by mlph \n";
        t0 = get_wall_time();
        MLPH mlph(method, b0, b1, g.size(), g.size(), g.degree_norm, adj_matrix, dual_value, g.size());
        mlph.run();
        duration = get_wall_time()-t0;
        cout << "time used: " << duration << "\n";
        time_duration_pricing_heur += duration;
        cout << "minimum reduced cost by mlph is: " << mlph.best_rc << endl;
        return mlph.best_rc;
    }

    bool CG::solve_mwis_gurobi(double cutoff, double& min_rc) {

        if (cutoff <= 0)
            return false;

        // setup the model now
        GRBEnv *env;
        vector<GRBVar> x;
        env = new GRBEnv();
        GRBModel model = GRBModel(*env);
        model.set(GRB_DoubleParam_TimeLimit, cutoff);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        model.set(GRB_StringAttr_ModelName, "MIP_MWIP");
        // Create variables and set them to be binary
        x.resize(g.size());
        for (int i = 0; i < g.size(); ++i){
            x[i] = model.addVar(0,1,0,GRB_BINARY);
        }
        model.update();

        // adjacent vertices cannot be selected simultaneously.
        for (int i = 0; i < g.size(); ++i){
            for (int j = i+1; j < g.size(); ++j){
                if (adj_matrix[i][j] == 1){
                    model.addConstr(x[i] + x[j] <= 1, "");
                }
            }
        }
        model.update();
        // cout << "max thread: " << thread_limit << "\n";
        model.set(GRB_IntParam_Threads, thread_limit);
        model.set(GRB_IntParam_Presolve, 0);
        // the objective
        GRBLinExpr tot = 1;
        for(int i = 0; i < g.size(); ++i){
            tot -= x[i] * dual_value[i];
        }
        model.setObjective(tot,GRB_MINIMIZE);
        model.update();
        auto t0 = get_wall_time();
        model.optimize();

        min_rc = model.get(GRB_DoubleAttr_ObjVal);
        bool mwis_optimal = (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL);
        
        if (mwis_optimal){
            optimal_mis.clear();
            for (int i = 0; i < g.size(); ++i){
                if (x[i].get(GRB_DoubleAttr_X) > 0.5){
                    optimal_mis.push_back(i);
                }
            }
            
            std::cout << "solved to optimality! time used: " << get_wall_time() - t0 << "\n";
        }else{
            std::cout << "reach cutoff time!\n";
        }

        delete env;
        return mwis_optimal;
    }

    void CG::test(int method, int column_selection, std::ofstream* output_file_sampling_stats, std::ofstream* output_file_cg_stats){
        
        double start_time = get_wall_time();
        int upper_col_limit = g.size();
        initializing_mis();
        min_reduced_cost_exact = -1.0;
        cg_iters = 0;

        double t0;
        double duration;
        bool mwis_optimal;
        while(true){
            cout << "iteration " << cg_iters++ << " : " << endl;
            t0 = get_wall_time();
            solve_restricted_master_problem(false);
            duration = get_wall_time()-t0;
            time_duration_master += duration;
            cout << "time used: " << duration << "\n";

            if (output_file_cg_stats!=nullptr){
                (*output_file_cg_stats) << cg_iters << "," 
                    << get_wall_time() - start_time << ","
                    << lp_bound << ",";
            }


            cout << "solve pricing mwis by heuristic \n";
            t0 = get_wall_time();
            Pricer* pricer = nullptr;
            
            if (method >= 3 && method <= 6){
                pricer = new MLPH(method, pricer_cutoff, g.size(), g.size(), g.degree_norm, adj_matrix, dual_value, upper_col_limit);            
            }else if (method == 7){
                pricer = new ACO(method, pricer_cutoff, g.size(), g.get_nb_edges(), 30, 
                                    g.degree_norm, adj_matrix, dual_value, upper_col_limit);
            }else if (method >= 8 && method <= 13){
                auto n_edges = g.get_nb_edges();
                pricer = new MWCP_solver(method, pricer_cutoff, dual_value, adj_list, adj_matrix, n_edges, upper_col_limit);
            }else{
                cout << "ERROR: at CG column_generation_heur()\n";
                assert(false);
            }
            pricer->column_selection = column_selection;
            pricer->run();                  
            if (output_file_sampling_stats!=nullptr){
                for (auto i = 0; i < pricer->niterations; i++){
                    (*output_file_sampling_stats) << cg_iters << "," << i+1 << "," 
                                    << pricer->num_neg_rc_current_iteration[i] << ","
                                    << pricer->best_rc_current_iteration[i] << "\n";
                }
            } 
            
            duration = get_wall_time()-t0;
            cout << "time used: " << duration << "\n";
            time_duration_pricing_heur += duration;
            cout << "minimum reduced cost by pricer is: " << pricer->best_rc << endl;

            if (pricer->num_neg_rc_col > 0){
                num_heur_runs_success++;
                cout << "# columns with negative reduced cost found by heuristic pricer: " << pricer->num_neg_rc_col << "\n";
                pricer->include_new_cols(mis_set, lp_vbasis);
                num_mis = mis_set.size();
            }
            else {
                cout << "heuristic pricer does not find any column with negative reduced cost!\n";
                t0 = get_wall_time();
                cout << "solve pricing mwis by exact tsm: \n";
                cout << "time budget used so far: " << (t0 - start_time) << "\n";
                auto mwis_cutoff = cutoff - (t0 - start_time);
                cout << "time budget remaining as cutoff time for solving mwis: " << cutoff - (t0 - start_time) << "\n";
                mwis_optimal = solve_mwis_tsm(mwis_cutoff, min_reduced_cost_exact);
                // mwis_optimal = solve_mwis_gurobi(1e8, min_reduced_cost_exact);
                cout << "minimum reduced cost by exact solver is: " << min_reduced_cost_exact << endl;
                duration = get_wall_time()-t0;
                cout << "time used: " << duration << "\n";
                time_duration_pricing_exact += duration;

                // add new columns
                mis_set.push_back(vector<int> (optimal_mis));
                num_mis++;
            }

            if (output_file_cg_stats!=nullptr){
                (*output_file_cg_stats) << pricer->num_neg_rc_col << ","
                            << pricer->best_rc << ","
                            << pricer->mean_rc << "," 
                            << pricer->median_rc << "," 
                            << pricer->stdev_rc << ","
                            << pricer->column_selection_time << ","
                            << pricer->dominate_selection_time << "\n";
            }
            
            if (pricer!=nullptr) delete pricer;
            if (!mwis_optimal) break;
            if (min_reduced_cost_exact > -0.000001) break;
            if (get_wall_time()-start_time > cutoff) break;
        }
        
        lp_optimal = mwis_optimal && min_reduced_cost_exact > -0.000001;
    }
}
