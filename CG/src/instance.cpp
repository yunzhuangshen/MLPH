#include "instance.h"
#include <cmath>
#include <random>
#include <iterator>
#include <algorithm>
#include "CG.h"
namespace GCP {
    static bool is_file_exist(const char *fileName)
    {
        std::ifstream infile(fileName);
        return infile.good();
    };

    Instance::Instance(string file_name, string input_dir, bool solve) : file_name{file_name}, input_dir{input_dir} {
        cout << file_name << endl;
        read_graph();

        string opt_file_name = input_dir + file_name + ".sol";
        if (solve){
            if (!is_file_exist(opt_file_name.c_str()))
                collect_train_data(opt_file_name);
        }
        if (is_file_exist(opt_file_name.c_str()))
            read_mis_sols(opt_file_name);
    }


    void Instance::read_mis_sols(std::string read_from){

        ifstream opt_file(read_from);
        string line;
        bool sol_val;
        double dual;
        int ctr = 0;

        std::vector<double> obj_coef;
        std::vector<bool> solution;
        while(!opt_file.eof()) {
            getline(opt_file, line);
            if (line == "EOF" || line == "-1" || line.size()==0) 
                break;

            std::stringstream stream(line);            
            if (ctr % 2 == 0 ) {  // read current dual values
                while (stream >> dual) 
                    obj_coef.push_back(dual);
                mis_obj_coefs.push_back(obj_coef);
            }else{ // read optimal solution
                while (stream >> sol_val) 
                    solution.push_back(sol_val);
                mis_sols.push_back(solution);
            }
            ctr++;
        }
        opt_file.close();
        ninst=ctr/2.;
    }

    void Instance::collect_train_data(std::string save_to){
        CG cg(*this, 1000, 1000, 1, 1314);

        vector<vector<double>> obj_coef;
        vector<std::vector<bool>> solution;
        cg.collect_training_data(obj_coef, solution);
        int num_mis_inst = obj_coef.size();

        ofstream opt_file(save_to);
        for (auto i = 0; i < num_mis_inst; i++){
            for (auto j = 0; j < n_nodes; j++){ 
                opt_file << obj_coef[i][j];
                if (j!=n_nodes-1) opt_file << " ";
            }
            opt_file << "\n";
            for (auto j = 0; j < n_nodes; j++){ 
                int sol_val = solution[i][j];
                opt_file << sol_val;
                if (j!=n_nodes-1) opt_file << " ";
            }
            opt_file << "\n";
        }
        opt_file.close();
    }

    void Instance::read_all_optimal_solution(){
        string opt_file_name = input_dir + file_name + ".allsol";
        ifstream opt_file(opt_file_name);
        if (!opt_file){
            cout << "optimal solution is not provided" << endl;
            return;
        }
        string line;
        bool READ_Point = 0;
        int node;
        int n0 = 0;
        int n1 = 1;
        while(!opt_file.eof()) {
            getline(opt_file, line);
            if (line == "EOF" || line == "-1") {
                break;
            }
            if (READ_Point){
                std::stringstream stream(line);
                while (stream >> node) {
                    optimal_value.push_back(node);
                    if (node == 0){
                        n0++;
                    }
                    if (node == 1){
                        n1++;
                    }
                }
            }
            if (line == "SET_INDEX"){
                READ_Point = 1;
            }
        }
        optimal_ratio = (float) n1 / (n0+n1);
        cout << "n0 is " << n0 << ", n1 is " << n1 << endl;

    }

    void Instance::read_optimal_solution(){
        string opt_file_name = input_dir + file_name + ".sol";
        ifstream opt_file(opt_file_name);
        if (!opt_file){
            cout << "optimal solution is not provided" << endl;
            return;
        }
        string line;
        vector<int> opt_sets;
        bool READ_Point = 0;
        int node;
        while(!opt_file.eof()) {
            getline(opt_file, line);
            if (line == "EOF" || line == "-1") {
                break;
            }
            if (READ_Point){
                std::stringstream stream(line);
                while (stream >> node) {
                    opt_sets.push_back(node);
                }
            }
            if (line == "SET_INDEX"){
                READ_Point = 1;
            }
        }
        int opt_obj = opt_sets.size();
        optimal_value = vector<bool>(n_sets, 0);
        int s;
        for (int i = 0; i < opt_sets.size(); ++i){
            s = opt_sets[i];
            optimal_value[s] = 1;
        }

        // verify if the solution is optimal
        vector<bool> covered_vertex(n_nodes, 0);
        int num_vertex_covered = 0;
        int v;
        for (int i = 0; i < opt_sets.size(); ++i){
            s = opt_sets[i];
            for (int j = 0; j < mis_sets[s].size(); ++j){
                v = mis_sets[s][j];
                if (covered_vertex[v] == 0){
                    covered_vertex[v] = 1;
                    num_vertex_covered++;
                }
            }
        }
        if (num_vertex_covered == n_nodes){
            cout << "optimal objective value is " << opt_obj << endl;
        } else{
            cout << "the mis sets provided are not feasible" << endl;
            assert(num_vertex_covered == n_nodes);
        }
    }

    void Instance::read_graph(){
        string input_file = input_dir  + file_name + ".col";
        ifstream file(input_file);
        string line, s1, s2;
        int v1, v2, ne = 0;
        int idx;

        while(!file.eof()) {
            getline(file, line);
//            cout << line << "\n";
            if (line[0] == 'p') {
                stringstream stream(line);
                stream >> s1 >> s2 >> n_nodes >> n_edges;
                adj_list.resize(n_nodes);
                cout << "number of nodes is " << n_nodes << "; number of edges is " << n_edges << endl;
            }
            if (line[0] == 'e'){
                stringstream stream(line);
                stream >> s1 >> v1 >> v2;
                if (v1 == v2) continue;
                if (find(adj_list[v1-1].begin(), adj_list[v1-1].end(), v2-1) == adj_list[v1-1].end()){
                    adj_list[v1-1].push_back(v2-1);
                    adj_list[v2-1].push_back(v1-1);
                    ne++;
                }
                // for (idx = 0; idx < adj_list[v1-1].size(); ++idx){
                //     if (adj_list[v1-1][idx] == v2 - 1){
                //         break;
                //     }
                // }
                // if (idx == adj_list[v1-1].size()){
                //     adj_list[v1-1].push_back(v2-1);
                //     adj_list[v2-1].push_back(v1-1);
                //     ne++;
                // }
            }
        }

        degree_norm = vector<double>(n_nodes);
        degree = vector<int>(n_nodes);

        max_node_degree_norm = 0.0;
        for (int i = 0; i < n_nodes; ++i){
            degree_norm[i] = (double) adj_list[i].size()/n_nodes;
            degree[i] = adj_list[i].size();
            if (max_node_degree_norm < degree_norm[i]){
                max_node_degree_norm = degree_norm[i];
            }
        }
//        cout << "Computing degree_norm done" << endl;
        n_edges = ne;
        cout << "number of undirected edges is " << n_edges << endl;
    }


    void Instance::compute_problem_specific_features(){
        set_size = vector<int>(n_sets, 0.0);
        max_set_degree_norm = vector<float>(n_sets, 0.0);
        ave_set_degree_norm = vector<float>(n_sets, 0.0);
        min_set_degree_norm = vector<float>(n_sets, 1.0);
        std_set_degree_norm = vector<float>(n_sets, 0.0);
        max_set_size = 0;
        int v;
        for (int i = 0; i < n_sets; ++i){
            set_size[i] = mis_sets[i].size();
            if (max_set_size < set_size[i]){
                max_set_size = set_size[i];
            }
            for (int j = 0; j < set_size[i]; ++j){
                v = mis_sets[i][j];
                ave_set_degree_norm[i] += degree_norm[v] / set_size[i];
                if (max_set_degree_norm[i] < degree_norm[v]){
                    max_set_degree_norm[i] = degree_norm[v];
                }
                if (min_set_degree_norm[i] > degree_norm[v]){
                    min_set_degree_norm[i] = degree_norm[v];
                }
            }
            for (int j = 0; j < set_size[i]; ++j){
                v = mis_sets[i][j];
                std_set_degree_norm[i] += pow(degree_norm[v] - ave_set_degree_norm[i], 2) / set_size[i];
            }
            std_set_degree_norm[i] = sqrt(std_set_degree_norm[i]);
        }


        vector<double> node_max_set_size(n_nodes, 0.0);
        vector<double> node_min_set_size(n_nodes, max_set_size);
        vector<double> node_ave_set_size(n_nodes, 0.0);
        int s;
        for (int i = 0; i < n_nodes; ++i){
            for (int j = 0; j < node_to_sets[i].size(); ++j){
                s = node_to_sets[i][j];
                node_ave_set_size[i] += (double) mis_sets[s].size() / node_to_sets[i].size();
                if (node_max_set_size[i] < mis_sets[s].size()){
                    node_max_set_size[i] = mis_sets[s].size();
                }
                if (node_min_set_size[i] > mis_sets[s].size()){
                    node_min_set_size[i] = mis_sets[s].size();
                }
            }
        }


        max_rel_set_size = vector<float>(n_sets, 0.0);
        ave_rel_set_size = vector<float>(n_sets, 0.0);
        min_rel_set_size = vector<float>(n_sets, 1.0);
        std_rel_set_size = vector<float>(n_sets, 0.0);
        double relative_size;
        for (int i = 0; i < n_sets; ++i){
            for (int j = 0; j < set_size[i]; ++j){
                v = mis_sets[i][j];
                relative_size = (double) set_size[i] / node_max_set_size[v];
                ave_rel_set_size[i] += relative_size / set_size[i];
                if (max_rel_set_size[i] < relative_size){
                    max_rel_set_size[i] = relative_size;
                }
                if (min_rel_set_size[i] > relative_size){
                    min_rel_set_size[i] = relative_size;
                }
            }
            for (int j = 0; j < set_size[i]; ++j){
                v = mis_sets[i][j];
                relative_size = (double) set_size[i] / node_max_set_size[v];
                std_rel_set_size[i] += pow(relative_size - ave_rel_set_size[i], 2) / set_size[i];
            }
            std_rel_set_size[i] = sqrt(std_rel_set_size[i]);
        }
    }

}