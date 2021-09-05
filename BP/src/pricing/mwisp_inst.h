
#ifndef __GCONVERT__
#define __GCONVERT__

#include <string>
#include <vector>
#include "../BP/probdata_coloring.h"
#include "../cg.h"

using namespace std;
class MWISP_INST{

public:

    METHOD_TYPE method_type;
    int nb_node;
    int nb_edge;
    vector<double> dual_values;
    vector<double> degree_norm;
    vector<vector<bool>> adj_matrix;
    vector<vector<int>> adj_list;

    MWISP_INST(TCLIQUE_GRAPH* graph, double* pi, METHOD_TYPE _method_type);


};
#endif
