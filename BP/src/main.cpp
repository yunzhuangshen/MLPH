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

/**@file   Coloring/src/main.c
 * @brief  Main file for C compilation
 * @author Gerald Gamrath
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/scipshell.h>
#include <string>
#include <vector>
#include <iostream>
#include <string>     // std::string, std::to_string
#include <boost/filesystem.hpp>
#include <iomanip>      // std::setprecision
#include "BP/reader_col.h"
#include "BP/branch_coloring.h"
#include "BP/coloringplugins.h"
#include "./BP/heur_init.h"
#include "cg.h"
#include "BP/heur_init.h"
using namespace std;

/**

command to run: ./BP [graph name] [heuristic pricer] [exact pricer] [branch-and-price cutoff time] [sample number]

heuristic pricer:
- 0 GREEDY
- 1 MSSP
- 2 SSSP
- 3 ACO
- 4 GUROBI
- 5 GUROBI_HEUR
- 6 TSM
- 7 FASTWCLQ
- 8 LSCC
- 9 None
- 10 MSSP_MIX
exact pricer:
- 0 tclique
- 1 tsm
- 2 gurobi
*/

string to_string_remove_trailing_zeros(double d){
  string s = to_string(d);
   s.erase(s.find_last_not_of('0')+1, std::string::npos);
   return s;
}

int
main(
   int                        argc,
   char**                     argv
   )
{
   string prob_dir = "../../GCB/";
   double cutoff_bp = 8000;
   double cutoff_pricing=30;
   EXACT_PRICER exact_pricer_type = static_cast<EXACT_PRICER>(0);

   string prob_name = argv[1]; // 1-FullIns_3
   string prob_path = prob_dir + prob_name + ".col";
   METHOD_TYPE method_type = static_cast<METHOD_TYPE>(stoi(argv[2]));
   int seed=stoi(argv[3]);

   double sample_factor = 10;
   double root_column_limit_factor = 1;
   double child_column_limit_factor = 0.1;
   string out_dir;
   if (method_type==METHOD_TYPE::BP_DEF){
      out_dir= "../results/BP_def/" + to_string(seed) + "/";
   } else if (method_type==METHOD_TYPE::BP_None){
      out_dir= "../results/BP_none/" + to_string(seed) + "/";
   } else if (method_type==METHOD_TYPE::BP_MLPH_FORCE_EXACT){
      out_dir= "../results/BP_MLPH_force_exact/" + to_string(seed) + "/";
   } else if (method_type==METHOD_TYPE::BP_MLPH){
      sample_factor = stod(argv[4]);
      root_column_limit_factor = stod(argv[5]);
      child_column_limit_factor = stod(argv[6]);
      out_dir= "../results/BP_MLPH_" + to_string_remove_trailing_zeros(sample_factor) + "_" 
         + to_string_remove_trailing_zeros(root_column_limit_factor) + "_" + to_string_remove_trailing_zeros(child_column_limit_factor) +"/" + to_string(seed) + "/";
   } else {
      cout << "method not recognized!\n"; assert(false);
   }
   boost::filesystem::create_directories(out_dir);

   SCIP_RETCODE retcode;
   SCIP* scip = NULL;
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeColoringPlugins(scip));
   SCIP_CALL( SCIPincludeHeurInit(scip, seed));
   SCIP_CALL( SCIPincludePricerColoring(scip, method_type, EXACT_PRICER::TCLIQUE, cutoff_bp, cutoff_pricing, sample_factor, root_column_limit_factor, child_column_limit_factor));
   retcode = SCIPreadProb(scip, prob_path.c_str(), NULL);

   switch( retcode )
   {
   case SCIP_NOFILE: SCIPinfoMessage(scip, NULL, "file <%s> not found\n", prob_path.c_str()); return SCIP_OKAY;
   case SCIP_PLUGINNOTFOUND: SCIPinfoMessage(scip, NULL, "no reader for input file <%s> available\n", prob_path.c_str()); return SCIP_OKAY;
   case SCIP_READERROR: SCIPinfoMessage(scip, NULL, "error reading file <%s>\n", prob_path.c_str()); return SCIP_OKAY;
   default: SCIP_CALL( retcode );
   } 

   // SCIPsetMessagehdlrQuiet(scip, true);
   SCIPsetIntParam(scip, "presolving/maxrestarts", 0);
   SCIPsetIntParam(scip, "randomization/permutationseed", seed);
   SCIPsetIntParam(scip, "randomization/randomseedshift", seed);
   SCIPsetRealParam(scip, "limits/time", cutoff_bp);
   SCIPsetIntParam(scip, "pricing/maxvars", 1e8);
   SCIPsetIntParam(scip, "pricing/maxvarsroot", 1e8);

   SCIP_CALL( SCIPsolve(scip) );

   // write gap curve
   ofstream output_file_gap_stats(out_dir+prob_name+".gap");
   for (auto i = 0; i < gaps.size(); i++){
      output_file_gap_stats << setprecision(5) << gap_times[i] << " " << gaps[i] << "\n";
   }
   output_file_gap_stats.close();


   // write solving statistics
   ofstream output_file_solving_stats(out_dir+prob_name+".txt");
   long long num_nodes = SCIPgetNNodes(scip);
   double primal_bound = SCIPgetPrimalbound(scip);
   double dual_bound = SCIPgetDualbound(scip);
   double opt_gap = SCIPgetGap(scip);
   double running_time = SCIPgetSolvingTime(scip);
   int solved = 0;
   SCIP_STATUS _status = SCIPgetStatus(scip);
   solved = (_status == SCIP_STATUS_OPTIMAL);
   output_file_solving_stats << solved << "," 
      << primal_bound << "," << dual_bound << "," 
      << opt_gap << "," << num_nodes << ","
      << running_time << "," << root_solving_time << "," << tot_exact_pricing_time << ","
      << tot_heuristic_pricing_time << "," << fixing_col_time << "," << primal_heur_time << "\n";
   output_file_solving_stats.close();

   SCIP_CALL( SCIPfree(&scip) );
   BMScheckEmptyMemory();
  if( retcode != SCIP_OKAY )
  {
     SCIPprintError(retcode);
     return -1;
  }

  return 0;
}
