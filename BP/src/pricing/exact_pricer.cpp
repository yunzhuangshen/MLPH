#include "exact_pricer.h"

namespace GCP {
    using namespace std;

    /* defines for rounding for tclique */
    #define MAXDNOM                1000LL
    #define MINDELTA               1e-03
    #define MAXDELTA               1e-09
    #define MAXSCALE               1000.0
    #define THRESHOLD               -1e-12
/** checks, whether the given scalar scales the given value to an integral number with error in the given bounds */


static
SCIP_Bool isIntegralScalar(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = EPSFLOOR(sval, 0.0); /*lint !e835*/
   upval = EPSCEIL(sval, 0.0); /*lint !e835*/

   return (SCIPrelDiff(sval, downval) <= maxdelta || SCIPrelDiff(sval, upval) >= mindelta);
}


/** get integral number with error in the bounds which corresponds to given value scaled by a given scalar;
 *  should be used in connection with isIntegralScalar()
 */
static
SCIP_Longint getIntegralVal(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;
   SCIP_Longint intval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = EPSFLOOR(sval, 0.0); /*lint !e835*/
   upval = EPSCEIL(sval, 0.0); /*lint !e835*/

   if( SCIPrelDiff(sval, upval) >= mindelta )
      intval = (SCIP_Longint) upval;
   else
      intval = (SCIP_Longint) downval;

   return intval;
}



/** generates improving variables using a stable set found by the algorithm for maximum weight clique,
 *  decides whether to stop generating cliques with the algorithm for maximum weight clique
 */

TCLIQUE_NEWSOL(tcliqueNewsolPricer)
{
    SCIP_PRICERDATA* pricerdata;
    int i;

    assert(acceptsol != NULL);
    assert(stopsolving != NULL);

    pricerdata = (SCIP_PRICERDATA*)tcliquedata;

    assert(pricerdata != NULL);
    assert(pricerdata->scip != NULL);
    assert(pricerdata->scalefactor > 0);

    *acceptsol = FALSE;
    *stopsolving = FALSE;



    /* if the stable set was already created in a former pricing round, we don't have to add it a second time */
    if ( !COLORprobStableSetIsNew(pricerdata->scip, cliquenodes, ncliquenodes) )
        return;

    vector<int> v(ncliquenodes);
    double obj = 0;
    for ( i = 0; i < ncliquenodes; i++ ){
        v[i] = cliquenodes[i];
        obj+= pricerdata->pi[cliquenodes[i]];
    }

        
    double rc = 1 - obj;
    if (rc < 0){
        pricerdata->improving_mwiss.push_back(v);
        if ( !pricerdata->solve_to_optimality ){
            *stopsolving = TRUE;
        }
    }
   /* accept the solution as new incumbent */
   *acceptsol = TRUE;
}/*lint !e715*/


    Exact_pricer::Exact_pricer(EXACT_PRICER _exact_type, double cutoff, MWISP_INST& inst, SCIP_PRICERDATA* pricerdata, double best_pricing_obj):
    mis_adj_list{inst.adj_list}, mis_adj_matrix{inst.adj_matrix}, cutoff{cutoff}, inst{inst}, pricerdata(pricerdata), best_pricing_obj{best_pricing_obj}
    {
        dual_values = inst.dual_values;
        exact_type = _exact_type;
        nb_node = dual_values.size();
        max_mwis_obj = -1e8;
    }

    void Exact_pricer::run(){
        switch (exact_type)
        {
        case EXACT_PRICER::EXACT_TSM: solve_tsm(); break;
        case EXACT_PRICER::EXACT_GUROBI: solve_gurobi(); break;
        case EXACT_PRICER::TCLIQUE: solve_tclique(); break;
        default:
            cout << "not recognizied option for exact solver: " << exact_type << "\n";
            exit(-1);
        }
    }

    void Exact_pricer::solve_tsm(){

        int nb_edge;
        long long** AdjacentList;
        long long* Node_Degree;
        long long* Node_Weight;
        long long* Node_Bound;
        
        nb_edge = (nb_node*(nb_node-1))/2. - inst.nb_edge;
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

        // solve the reduced graph using TSM
        TSM_SOLVER::TSM solver;
        char const *name2 = "ret";
        solver.tsm("", name2, nb_node, nb_edge, cutoff, AdjacentList, Node_Degree, Node_Weight, Node_Bound);
        objs = solver.sol_objs; sols = solver.sols;
        isOptimal = solver.isOptimal;
        free(Node_Degree);
        free(Node_Bound);
        free(Node_Weight);
        for (auto i = 0; i < nb_node; i++)
            free(AdjacentList[i]);
        free(AdjacentList);

        for (auto i = 0; i < sols.size(); i++){
            double obj = 0;
            double rc = 1;
            for (auto v : sols[i]){
                obj += dual_values[v];
            }

            if (obj > max_mwis_obj)
                max_mwis_obj = obj;
            rc = 1 - obj;
            if (rc < THRESHOLD){
                neg_rc_vals.push_back(rc);
                neg_rc_cols.push_back(sols[i]);
            }
            
        }
    }

    void Exact_pricer::solve_tclique(){
        TCLIQUE_GRAPH*   cgraph;                /* the complementary graph, used for tclique-algorithm */
        bool        weightsIntegral;
        unsigned int        scalesuccess;
        int*             maxstablesetnodes;     /* pointer to store nodes of the maximum weight clique */
        int              nmaxstablesetnodes;    /* number of nodes in the maximum weight clique */
        TCLIQUE_WEIGHT   maxstablesetweight;    /* weight of the maximum weight clique */
        TCLIQUE_STATUS   status;                /* status of clique-computation */
        pricerdata->improving_mwiss.clear();
        /* get the complementary graph from the current cons */
        cgraph = COLORconsGetComplementaryGraph(pricerdata->scip);
        SCIPallocBufferArray(pricerdata->scip, &maxstablesetnodes, nb_node);
        /* get dual solutions and set weight of nodes */
        weightsIntegral = TRUE;
        for (int i = 0; i < nb_node; i++ )
        {
            pricerdata->pi[i] = SCIPgetDualsolSetppc(pricerdata->scip, pricerdata->constraints[i]);

            if( !isIntegralScalar(pricerdata->pi[i], 1.0, -MINDELTA, MAXDELTA) )
            {
                weightsIntegral = FALSE;
            }
        }
        /* are weigths integral? */
        if( weightsIntegral )
        {
            pricerdata->scalefactor = 1.0;
            scalesuccess = TRUE;
        }
        else
        {
            /* compute factor, which makes the weights integral */
            scalesuccess = FALSE;
            SCIPcalcIntegralScalar(pricerdata->pi, nb_node, -MINDELTA, MAXDELTA, MAXDNOM, MAXSCALE,
                &(pricerdata->scalefactor), &scalesuccess);

        }
        // std::cout << "scale factor: " << pricerdata->scalefactor << "\n";
        assert(scalesuccess);
        /* change the weights for the nodes in the graph to the dual solution value * scalefactor */
        for (int i = 0; i < nb_node; i++ )
            tcliqueChangeWeight(cgraph, i, getIntegralVal(pricerdata->pi[i], pricerdata->scalefactor, -MINDELTA, MAXDELTA)); /*lint !e712 !e747*/
        int force_mlhp = best_pricing_obj != 1;
        best_pricing_obj =  getIntegralVal(best_pricing_obj, pricerdata->scalefactor, -MINDELTA, MAXDELTA);
        /* compute maximal clique */
        tcliqueMaxClique(NULL, NULL, NULL, NULL, cgraph, tcliqueNewsolPricer, (TCLIQUE_DATA*)pricerdata, maxstablesetnodes,
            &(nmaxstablesetnodes), &maxstablesetweight, 0,
            best_pricing_obj, pricerdata->maxtcliquenodes, 0, INT_MAX, -1,
            NULL, &status, cutoff, force_mlhp);
        SCIPfreeBufferArray(pricerdata->scip, &maxstablesetnodes);
        
        // if (pricerdata->solve_to_optimality) assert(status == TCLIQUE_OPTIMAL);
        // else assert(status == TCLIQUE_USERABORT);
        isOptimal=status==TCLIQUE_OPTIMAL;

        for (int i = 0; i < pricerdata->improving_mwiss.size(); i++){
            double obj=0; 
            double rc = 1;
            for (auto v : pricerdata->improving_mwiss[i]){
                obj += dual_values[v];
            }

            if (obj > max_mwis_obj)
                max_mwis_obj = obj;
            rc = 1 - obj;
            neg_rc_vals.push_back(rc);
            neg_rc_cols.push_back(pricerdata->improving_mwiss[i]);
        }
        pricerdata->improving_mwiss.clear();
    }


    void Exact_pricer::solve_gurobi() {

        // setup the model now
        GRBEnv *env;
        vector<GRBVar> x;
        env = new GRBEnv();
        GRBModel model = GRBModel(*env);
        model.set(GRB_IntParam_PoolSolutions, 1e8); //retain all solutions
        model.set(GRB_DoubleParam_TimeLimit, 1e8);
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

        // the objective
        GRBLinExpr tot = 0;
        for(int i = 0; i < nb_node; ++i){
            tot += x[i] * dual_values[i];
        }
        model.setObjective(tot,GRB_MAXIMIZE);
        model.update();
        model.optimize();
        isOptimal = model.get(GRB_IntAttr_Status)==GRB_OPTIMAL;
        auto nsol = model.get(GRB_IntAttr_SolCount);
        for (int i = 0; i < nsol; i++){
            vector<int> sol;
            model.set(GRB_IntParam_SolutionNumber, i);
            auto obj = model.get(GRB_DoubleAttr_PoolObjVal);
            if (obj > max_mwis_obj)
                max_mwis_obj = obj;

            auto rc = 1. - obj;
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
}
