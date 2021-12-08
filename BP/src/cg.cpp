/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pricer_coloring.c
 * @brief  variable pricer for the vertex coloring problem
 * @author Gerald Gamrath
 *
 * This file implements the pricer for the coloring algorithm.
 *
 * It computes maximal stable sets in the current graph whose corresponding variables can improve
 * the current LP solution.  This is done by computing a maximum weighted stable set in the current
 * graph with dual-variables of the node constraints as weights. A variable can improve the
 * solution, if the weight of the corresponding stable set is larger than 1, since it then has
 * negative reduced costs, which are given by (1 - weight of the set).
 *
 * The pricer first tries to compute such a stable set using a a greedy-method. If it fails, the tclique-algorithm is
 * used on the complementary graph. This is a branch-and-bound based algorithm for maximal cliques,
 * included in SCIP.  In this case, not only the best solution is added to the LP, but also all other
 * stable sets found during the branch-and-bound process that could improve the current LP solution
 * are added, limited to a maximal number that can be changed by a parameter.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "cg.h"
#include "BP/reader_col.h"
#include "BP/cons_storeGraph.h"
#include "pricing/MLPH.h"
#include "pricing/heur_pricer_misc.h"
#include "pricing/exact_pricer.h"
#include "pricing/mwisp_inst.h"
#include <cmath>

#define PRICER_NAME            "coloring"
#define PRICER_DESC            "pricer for coloring"
#define PRICER_PRIORITY        5000000
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */


/* default values for parameters */
#define DEFAULT_OPTIMALITY        TRUE
#define DEFAULT_MAXROUNDSROOT   -1
#define DEFAULT_MAXROUNDSNODE   -1
#define DEFAULT_MAXTCLIQUENODES INT_MAX

using namespace std;

double tot_exact_pricing_time = 0;
double tot_heuristic_pricing_time = 0;
double root_solving_time = 0;
double root_time_start = 0;
int enter_root=0;
vector<double> gaps;
vector<double> gap_times;

/*
 * Local methods
 */

/** returns whether the graph has an uncolored node
 */
static
SCIP_Bool hasUncoloredNode(
   TCLIQUE_GRAPH*        graph,              /**< the graph that should be colored */
   SCIP_Bool*            colored             /**< array of booleans, colored[i] == TRUE iff node i is colored */
   )
{
   int i;

   assert(graph != NULL);
   assert(colored != NULL);

   for ( i = 0; i < tcliqueGetNNodes(graph); i++)
   {
      /* node not yet colored */
      if (!colored[i])
      {
	return TRUE;
      }
   }
   return FALSE;
}

static double get_wall_time(){
   struct timeval time;
   if (gettimeofday(&time,NULL)){
      return 0;
   }
   return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

/** For farkas pricing, computes a stable set with a greedy-method.  attention: the weight of the maximum stable set is not computed! */
static
SCIP_RETCODE greedyStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   TCLIQUE_GRAPH*        graph,              /**< pointer to graph data structure */
   SCIP_Bool*            colored,            /**< array for marking yet colored nodes */
   int*                  maxstablesetnodes,  /**< pointer to store nodes of the maximum weight stableset */
   int*                  nmaxstablesetnodes  /**< pointer to store number of nodes in the maximum weight stableset */
   )
{
   SCIP_Bool indnode;
   int nnodes;
   int i;
   int j;
   int* degrees;
   int* sortednodes;
   SCIP_Real* values;    /* values for sorting the nodes: deg(v)+w(v)*nnodes  */

   assert(scip != NULL);
   assert(graph != NULL);
   assert(maxstablesetnodes != NULL);
   assert(nmaxstablesetnodes != NULL);

   /* get number of nodes */
   nnodes = tcliqueGetNNodes(graph);
   *nmaxstablesetnodes = 0;

   /* get the  degrees for the nodes in the graph */
   degrees = tcliqueGetDegrees(graph);
   SCIP_CALL( SCIPallocBufferArray(scip, &values, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortednodes, nnodes) );

   /* set values to the nodes which are used for sorting them */
   /* value = degree of the node + weight of the node * number of nodes, therefore the yet colored nodes
      (which have weight 0) have lower values than the not yet colored nodes which have weight 1 */
   for ( i = 0; i < nnodes; i++ )
   {
      sortednodes[i] = i;
      values[i] = ( colored[i] == TRUE ? degrees[i] : degrees[i]+nnodes );
   }

   /* sort the nodes w.r.t. the computed values */
   SCIPsortDownRealInt(values, sortednodes, nnodes);

   /* insert first node */
   maxstablesetnodes[0] = sortednodes[0];
   (*nmaxstablesetnodes) = 1;
   for ( i = 1; i < nnodes; i++)
   {
      /* check whether node is independent to nodes in the set */
      indnode = TRUE;
      for ( j = 0; j < (*nmaxstablesetnodes); j++ )
      {
         if ( tcliqueIsEdge(graph, sortednodes[i], maxstablesetnodes[j]) )
         {
            indnode = FALSE;
            break;
         }
      }
      if ( indnode == TRUE )
      {
         /* node is independent, thus add it to the set */
         maxstablesetnodes[*nmaxstablesetnodes] = sortednodes[i];
         (*nmaxstablesetnodes) = (*nmaxstablesetnodes)+1;
      }

   }
   SCIPfreeBufferArray(scip, &sortednodes);
   SCIPfreeBufferArray(scip, &values);

   return SCIP_OKAY;
}


/*
 * Callback methods of variable pricer
 */

/** copy method for pricer plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRICERCOPY(pricerCopyColoring)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(pricer != NULL);
   assert(strcmp(SCIPpricerGetName(pricer), PRICER_NAME) == 0);

   return SCIP_OKAY;
}


/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreeColoring)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   /* get pricerdata */
   pricerdata = SCIPpricerGetData(pricer);

   /* free memory for pricerdata*/
   if ( pricerdata != NULL )
   {
      delete pricerdata;
      // SCIPfreeBlockMemory(scip, &pricerdata);
   }

   SCIPpricerSetData(pricer, NULL);
   return SCIP_OKAY;
}



/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
static
SCIP_DECL_PRICERINITSOL(pricerInitsolColoring)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   pricerdata->bbnode = NULL;

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(pricerdata->pi), COLORprobGetNNodes(scip)) );
   return SCIP_OKAY;
}



/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolColoring)
{
   SCIP_PRICERDATA* pricerdata;
   int i;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   SCIPfreeBlockMemoryArray(scip, &(pricerdata->pi), COLORprobGetNNodes(scip));

   return SCIP_OKAY;
}




/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostColoring)
{
   SCIP_PRICERDATA* pricerdata;            /* the data of the pricer */

   TCLIQUE_GRAPH*   graph;                 /* the current graph */
   int              nnodes;                /* number of nodes in the graph */
   SCIP_Real        mwis_obj;/* weigth of the maximal stable set computed by the greedy */
   /* variables used for scaling the rational dual-solutions to integers */
   SCIP_Bool        scalesuccess;

   int              i;
   int              j;

   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   // std::cout << COLORprobGetNStableSets(pricerdata->scip) << "\n";
   // std::cout << "current node: " << SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) << "\n";
   /* count down number of remaining pricing rounds at the current node */
   if (enter_root==0){
      root_time_start = SCIPgetSolvingTime(scip);
      gaps.clear();
      gap_times.clear();
      enter_root = 1;
   }else if (enter_root == 1 && 
      SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) != SCIPnodeGetNumber(SCIPgetRootNode(scip))){
      
      root_solving_time = SCIPgetSolvingTime(scip) - root_time_start;
      enter_root = 2;
   }

   if ( pricerdata->bbnode == SCIPgetCurrentNode(scip) )
   {
      if ( pricerdata->noderounds > 0 )
         pricerdata->noderounds--;
   }
   else
   {
      if ( pricerdata->bbnode == NULL )
      {
         pricerdata->noderounds = pricerdata->maxroundsroot;
         pricerdata->lowerbound = - SCIPinfinity(scip);
      }
      else
      {
         pricerdata->noderounds = pricerdata->maxroundsnode;
         pricerdata->lowerbound = - SCIPinfinity(scip);
      }
      pricerdata->bbnode = SCIPgetCurrentNode(scip);
   }
   /* stop pricing if limit for pricing rounds reached */
   if ( pricerdata->noderounds == 0 )
   {
      SCIPdebugMessage("maxrounds reached, pricing interrupted\n");

      /* set result and lowerbound pointer */
      *result = SCIP_DIDNOTRUN;
      *lowerbound = pricerdata->lowerbound;
      
      return SCIP_OKAY;
   }

   //record gap and time
   double lb = SCIPgetLowerbound(scip);
   double ub = SCIPgetUpperbound(scip);
   if (lb != -SCIPinfinity(scip)){
      lb = ceil(lb);
      double gap = ((ub - lb)/ub)*100;
      if (gaps.size()==0){
         gaps.push_back(gap);
         gap_times.push_back(SCIPgetSolvingTime(scip));
      }else if (gap < gaps[gaps.size()-1]){
         gaps.push_back(gap);
         gap_times.push_back(SCIPgetSolvingTime(scip));
      }
   }

   /* set result pointer */
   *result = SCIP_SUCCESS;

   /* get graph and number of nodes */
   graph = COLORconsGetCurrentGraph(scip);
   assert(graph != NULL);
   nnodes = tcliqueGetNNodes(graph);
   // int nedges = tcliqueGetNEdges(graph);
   // cout << "stats: " << nnodes << " " << nedges << "\n";
   // SCIPinterruptSolve(scip);
   // *result = SCIP_DIDNOTRUN;
   // return SCIP_OKAY;
   assert(SCIPgetNVars(scip) == COLORprobGetNStableSets(scip));

   /* get constraints */
   pricerdata->constraints = COLORprobGetConstraints(scip);

   /* get dual solutions and save them in pi */
   for ( i = 0; i < nnodes; i++){
      pricerdata->pi[i] = SCIPgetDualsolSetppc(scip, pricerdata->constraints[i]);
   }
   MWISP_INST inst(graph, pricerdata->pi, pricerdata->method_type);

   /* ......heuristic........ */
   GCP::Pricer* my_pricer = nullptr;
   bool is_root = SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == SCIPnodeGetNumber(SCIPgetRootNode(scip));

   if (pricerdata->method_type==METHOD_TYPE::BP_MLPH || 
         pricerdata->method_type==METHOD_TYPE::BP_MLPH_FORCE_EXACT){
         my_pricer = new GCP::MLPH(is_root, pricerdata->method_type, pricerdata->cutoff_pricing_heur, inst, pricerdata->sample_factor);            
   }else if (
        pricerdata->method_type==METHOD_TYPE::BP_DEF||
        pricerdata->method_type==METHOD_TYPE::BP_None 
        ){
      my_pricer = new GCP::Heur_pricer_misc(pricerdata->method_type, pricerdata->cutoff_pricing_heur, inst);         
   }else{
         cout << "ERROR: at Reduce column_generation_heur()\n";
         assert(false);
   }

   // upon selecting a new node, reset the stats
   if (pricerdata->cur_node_id!=SCIPnodeGetNumber(SCIPgetCurrentNode(scip))){
      pricerdata->cur_node_id = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
      pricerdata->low_lp_obj_improve_count = 0;
   
   }else{ // contiuing solving current node
      double lp_obj_diff = pricerdata->lp_obj - SCIPgetLPObjval(scip);
      if (lp_obj_diff < 0.05) pricerdata->low_lp_obj_improve_count++;
      else pricerdata->low_lp_obj_improve_count = 0;
      // cout << "lp_obj_diff: " << lp_obj_diff << "\n";
      // cout << "low_lp_obj_improve_count: " << pricerdata->low_lp_obj_improve_count << "\n";
   }
   // record current lp objective value
   pricerdata->lp_obj = SCIPgetLPObjval(scip);

   double t0 = get_wall_time();
   my_pricer->run();
   double heur_time = get_wall_time() - t0;     
   tot_heuristic_pricing_time += heur_time;
   bool use_exact;

   // force SSSP to use exact to tackle tail off effect
   // cout << "LP objective: " << SCIPgetLPObjval(scip) << "\n";
   bool force_mlph = false;
   double best_pricing_obj = 1;
   vector<vector<int>> neg_rc_cols;
   vector<double> neg_rc_vals;

   force_mlph = pricerdata->method_type==METHOD_TYPE::BP_MLPH_FORCE_EXACT &&
       !my_pricer->neg_rc_cols.empty() &&
       pricerdata->low_lp_obj_improve_count >= 5 &&
       pricerdata->low_lp_obj_improve_count % 5 == 0; 
   
   if (force_mlph){
      best_pricing_obj = my_pricer->best_obj;
      neg_rc_cols = my_pricer->neg_rc_cols;
      neg_rc_vals = my_pricer->neg_rc_vals;
      cout << "force mlph!\n";
   }


   use_exact = force_mlph || my_pricer->neg_rc_cols.empty();

   t0 = get_wall_time();
   /* solve using an exact method */ 
   bool is_exact_optimal = false;  
   if (use_exact){
      delete my_pricer;
      double t = pricerdata->bp_cutoff - SCIPgetSolvingTime(scip);
      if (t<=0) t=0.001;
      // std::cout << "use_exact: " << t << "\n";

      my_pricer = new GCP::Exact_pricer(t, inst, pricerdata, best_pricing_obj);
      my_pricer->run();
      GCP::Exact_pricer* exact_pricer = dynamic_cast<GCP::Exact_pricer*>(my_pricer);
      is_exact_optimal = exact_pricer->isOptimal;
      if (exact_pricer->isOptimal && 
            SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL && 
            SCIPisFeasGT(scip, exact_pricer->max_mwis_obj, 1.0) )
      {
         pricerdata->lowerbound = MAX( pricerdata->lowerbound,
            SCIPgetLPObjval(scip)/exact_pricer->max_mwis_obj);
      }
   }

   double exact_time = get_wall_time() - t0;
   // - use lp lower bound
   if (ceil(pricerdata->lowerbound) == ceil(SCIPgetLPObjval(scip))){
      SCIPupdateLocalLowerbound(scip, ceil(SCIPgetLPObjval(scip)));
      *stopearly = TRUE;
      // cout << "early stopping by lagrangian bound\n";
   }

   tot_exact_pricing_time += get_wall_time() - t0;
   
   if (*stopearly == TRUE){
      *result = SCIP_DIDNOTRUN;
   }else if (my_pricer->neg_rc_cols.size()>0){
      if (force_mlph){
         my_pricer->neg_rc_cols.insert(my_pricer->neg_rc_cols.end(),neg_rc_cols.begin(), neg_rc_cols.end());
         my_pricer->neg_rc_vals.insert(my_pricer->neg_rc_vals.end(),neg_rc_vals.begin(), neg_rc_vals.end());
      }
      my_pricer->add_columns(pricerdata, is_root);

   }

   delete my_pricer;

   return SCIP_OKAY;
}/*lint !e715*/


/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasColoring)
{
   TCLIQUE_GRAPH* graph;
   int nnodes;                  /* number of nodes */
   int* maxstablesetnodes;      /* array containig the nodes of the max stable set */
   int nmaxstablesetnodes;      /* number of nodes in stable set */
   int setnumber;               /* number of already found stable sets */
   SCIP_VAR* var;               /* var for the actual stable set */
   SCIP_CONS** constraints;     /* array of added constraints */
   SCIP_Bool* colored;          /* array for marking of yet colored nodes
                                   colored_i = true iff node i is already colored */
   int**              stablesets;
   int*               nstablesetelements;
   int                nstablesets;
   int i;
   int j;
   assert(scip != NULL);

   graph = COLORconsGetCurrentGraph(scip);
   assert(graph != NULL);

   nnodes = COLORprobGetNNodes(scip);
   assert(nnodes > 0);

   /* get the node-constraits */
   constraints = COLORprobGetConstraints(scip);
   assert(constraints != NULL);

   /* get all yet computed stable sets */
   COLORprobGetStableSets(scip, &stablesets, &nstablesetelements, &nstablesets);
   assert(stablesets != NULL && nstablesetelements != NULL);
   assert(nstablesets >= 0);
   assert(nnodes == tcliqueGetNNodes(graph));

   /* allocate memory for arrays */
   SCIP_CALL( SCIPallocBufferArray( scip, &colored, nnodes) );
   SCIP_CALL( SCIPallocBufferArray( scip, &maxstablesetnodes, nnodes) );
   nmaxstablesetnodes = 0;

   /* fill colored-array with FALSE */
   BMSclearMemoryArray(colored, nnodes);

   /* go through all stable sets and set colored to true for nodes in them */
   for ( i = 0; i < nstablesets; i++ )
   {
      if ( !SCIPisFeasZero(scip, SCIPvarGetUbLocal( COLORprobGetVarForStableSet(scip, i)))
         && (SCIPgetNNodes(scip) == 0 || SCIPvarIsInLP(COLORprobGetVarForStableSet(scip, i))
            || SCIPgetRootNode(scip) == SCIPgetCurrentNode(scip) ) )
      {
         for ( j = 0; j < nstablesetelements[i]; j++ )
         {
            colored[stablesets[i][j]] = TRUE;
         }
      }
   }

   /* create maximal Stable Sets until all Nodes are covered */
   while ( hasUncoloredNode(graph, colored) )
   {
      SCIP_CALL( greedyStableSet(scip, graph, colored, maxstablesetnodes, &nmaxstablesetnodes) );
      SCIPsortDownInt(maxstablesetnodes, nmaxstablesetnodes);
      SCIP_CALL( COLORprobAddNewStableSet(scip, maxstablesetnodes, nmaxstablesetnodes, &setnumber) );
      assert(setnumber != -1);

      /* create variable for the stable set and add it to SCIP */
      SCIP_CALL( SCIPcreateVar(scip, &var, NULL, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY,
            TRUE, TRUE, NULL, NULL, NULL, NULL, (SCIP_VARDATA*) (size_t) setnumber) ); /*lint !e571*/
      SCIP_CALL( COLORprobAddVarForStableSet(scip, setnumber, var) );
      SCIPvarMarkDeletable(var);
      SCIP_CALL( SCIPaddPricedVar(scip, var, 1.0) );
      SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );

      for ( i = 0; i < nmaxstablesetnodes; i++ )
      {
         /* add variable to node constraints of nodes in the set */
         SCIP_CALL( SCIPaddCoefSetppc(scip, constraints[maxstablesetnodes[i]], var) );
         /* mark node as colored */
         colored[maxstablesetnodes[i]] = TRUE;
      }
   }
   /* free memory */
   SCIPfreeBufferArray(scip, &maxstablesetnodes);
   SCIPfreeBufferArray(scip, &colored);

   return SCIP_OKAY;
}/*lint !e715*/


/*
 * variable pricer specific interface methods
 */

/** creates the coloring variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerColoring(
   SCIP*                 scip,                /**< SCIP data structure */
   METHOD_TYPE           method_type,
   double                cutoff_bp,
   double                cutoff_pricing_heur, 
   double                sample_factor, 
   double                root_column_limit_factor, 
   double                child_column_limit_factor
   )
{
   SCIP_PRICERDATA* pricerdata;
   SCIP_PRICER* pricer;
   pricerdata = new SCIP_PRICERDATA();
   // SCIP_CALL( SCIPallocBlockMemory(scip, &pricerdata) );
   pricerdata->scip = scip;
   pricerdata->method_type = method_type;
   pricerdata->sample_factor = sample_factor;
   pricerdata->root_column_limit_factor = root_column_limit_factor;
   pricerdata->child_column_limit_factor = child_column_limit_factor;

   pricerdata->cur_node_id = -1;
   pricerdata->lp_obj = -1;
   pricerdata->low_lp_obj_improve_count = 0;
   tot_exact_pricing_time = 0.;
   tot_heuristic_pricing_time = 0.;
   pricer = NULL;
   /* include variable pricer */
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostColoring, pricerFarkasColoring, pricerdata) );
   assert(pricer != NULL);

   /* include non-fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPricerCopy(scip, pricer, pricerCopyColoring) );
   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreeColoring) );
   SCIP_CALL( SCIPsetPricerInitsol(scip, pricer, pricerInitsolColoring) );
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolColoring) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "pricers/coloring/solve_to_optimality",
         "For the tclique, should the best variables be addded to the problem instead of adding the first found variables?",
         &pricerdata->solve_to_optimality, FALSE, DEFAULT_OPTIMALITY, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "pricers/coloring/maxroundsroot",
         "maximum number of pricing rounds in the root node (-1: no limit)",
         &pricerdata->maxroundsroot, TRUE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "pricers/coloring/maxroundsnode",
         "maximum number of pricing rounds in each node (except root node)(-1: no limit)",
         &pricerdata->maxroundsnode, TRUE, DEFAULT_MAXROUNDSNODE, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "pricers/coloring/maxtcliquenodes",
         "maximum number of B&B-nodes used in the tclique-algorithm",
         &pricerdata->maxtcliquenodes, TRUE, DEFAULT_MAXTCLIQUENODES, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "pricers/coloring/cutoff_pricing_heur",
         "the cutoff time of the heur_pricer",
         &pricerdata->cutoff_pricing_heur, TRUE, cutoff_pricing_heur, 0, 1.e8, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "pricers/coloring/bp_cutoff",
         "the cutoff time of the heur_pricer",
         &pricerdata->bp_cutoff, TRUE, cutoff_bp, 0, 1.e8, NULL, NULL) );

   return SCIP_OKAY;
}
