#include "mwisp_inst.h"


MWISP_INST::MWISP_INST(TCLIQUE_GRAPH* graph, double* pi, METHOD_TYPE _method_type){
    nb_node =  tcliqueGetNNodes(graph);
    nb_edge = tcliqueGetNEdges(graph);
    int* degrees = tcliqueGetDegrees(graph);
    dual_values.resize(nb_node);
    adj_list = vector<vector<int>>(nb_node, vector<int>());

    for (auto i = 0; i < nb_node; i++){
        dual_values[i]=pi[i];
        adj_list[i].resize(degrees[i]);
    }

    // if (pricer_type== HEUR_PRICER::SSSP || 
    //     pricer_type==HEUR_PRICER::MSSP || 
    //     pricer_type==HEUR_PRICER::ACO ||
    //     pricer_type==HEUR_PRICER::GUROBI || 
    //     pricer_type==HEUR_PRICER::GUROBI_HEUR||
    //     pricer_type==HEUR_PRICER::GREEDY)
    // {
        adj_matrix = vector<vector<bool>>(nb_node, vector<bool>(nb_node, false));

        for (auto i = 0; i < nb_node; i++){
            int ctr=0;
            for (auto j = 0; j < nb_node; j++){
                if (i!=j && tcliqueIsEdge(graph, i, j)){
                    adj_matrix[i][j]=true;
                    adj_list[i][ctr++]=j;
                }
            }
            assert(ctr==degrees[i]);
        }

        degree_norm.clear();
        degree_norm.assign(degrees, degrees+nb_node);
        int degree_max = 0;
        for (auto i = 0; i < nb_node; i++)
            if (degrees[i]>degree_max) degree_max=degree_norm[i];
        for (auto i = 0; i < nb_node; i++)
            degree_norm[i] = degrees[i]/degree_max;
    // } 
    // else if (pricer_type== HEUR_PRICER::TSM || 
    //     pricer_type==HEUR_PRICER::LSCC || 
    //     pricer_type==HEUR_PRICER::FASTWCLQ){
    // }

}