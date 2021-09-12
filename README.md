# A Machine-Learning-Based Pricing Heuristic for Column Generation

This repository holds the code for using brand-and-price to solve the graph coloring problem. We propose a machine-learning-based heuristic pricing method to accelarate the progress of column generation. The code is organized as follows:

- GCB, containing Graph Coloring Benchmarks
- CG, containing c++ code for column generation and python code for analyzing the results.
- BP, containing c++ code for branch-and-price and python code for analyzing the results.

## Requirements
The C++ code can then be built with cmake (version >= 3.10) with:
- openMP - https://www.openmp.org/
- Boost - https://www.boost.org/
- Gurobi (9.0.1) - https://www.gurobi.com/ (Acdemic License is Required)
- SCIP (7.0.2) - https://www.scipopt.org/

The python code requires:
- numpy 
- scipy 
- bayesian-optimization - https://github.com/fmfn/BayesianOptimization 

## Run scrips to reproduce results:
1. python3 01-train-and-optimize.py
2. python3 02-cg.py (nCPUs $\in [4,8,12...]$)
3. python3 03-bp.py (nCPUs $\in [1,2,3,...]$)


For the second and third step, you can specificy the number of available CPUs in the python script.

## Results
The results are in the two newly created folders:
- `results_cg' contains the results for column generation
- `results_bp' containing the results for branch-and-price

The Figures and Tables in our main paper corresonponds to the results files respectively: 
- data for Figure 2: 
    - 'results_cg/small/lp-curve'
    - 'results_cg/small/solving-curve'
- data for Figure 3:
    - 'results_cg/small/compare_figure.txt'
    - 'results_cg/small/compare_number.txt'
- data for Figure 4:
    - 'results_cg/cs-large/lp-curve-cg'
    - 'results_cg/cs-large/lp-cg'
- data for Figure 5:
    - 'results_bp/gap_curve_BP_MLPH_10._1._0.1-BP_def'
- Table 2:
    - 'results_cg/large/table_solving_stats.tex'
- Table 3:
    - 'results_cg/large/table_rc.tex'  
- Table 4-6:
    - 'results_bp/table_BP_MLPH_10._1._0.1-BP_def/time_for_all_solved/*.tex'
    - 'results_bp/table_BP_MLPH_10._1._0.1-BP_def/gap_for_all_not_solved/*.tex'
    - 'results_bp/table_BP_MLPH_10._1._0.1-BP_def/number_solve_for_not_all_solved/*.tex'
