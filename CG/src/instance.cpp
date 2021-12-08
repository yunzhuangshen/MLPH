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
            load_training_data(opt_file_name);
    }


    void Instance::load_training_data(std::string read_from){

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
        nb_pp=ctr/2.;
        opt_file.close();
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
        n_edges = ne;
        cout << "number of undirected edges is " << n_edges << endl;
    }

}