#ifndef INSTANCE_H
#define INSTANCE_H
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <random>
#include <iterator>
#include <assert.h>
using namespace std;

namespace GCP {
    class Instance {
        int n_nodes;
        int n_edges;
        const string file_name;
        const string input_dir;
        void read_graph();
    
    public:
        int nb_pp = 0;
        vector<vector<int>> adj_list;
        vector<double> degree_norm;
        float max_node_degree_norm;
        vector<int> degree;
        std::vector<std::vector<double>> mis_obj_coefs;
        std::vector<std::vector<bool>> mis_sols;
        void load_training_data(std::string read_from);
        void collect_train_data(std::string save_to);
        explicit Instance(string file_name, string input_dir, bool _solve);

        // Size of the graph (i.e. number of nodes).
        int size() const { return n_nodes; }
        int get_nb_edges() const { return n_edges; }
        string get_file_name() const { return file_name; }
        // Prints out the cost matrix of the graph.
        friend ostream& operator<<(ostream& out, const Instance& g);
    };
}

#endif
