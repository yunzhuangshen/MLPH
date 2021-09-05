# A Machine-Learning-Based Pricing Heuristic for Column Generation

This repository holds the code for using brand-and-price to solve the graph coloring problem. We propose a machine-learning-based heuristic pricing method to accelarate the progress of column generation. The code is organized as follows:

- GCB, containing Graph Coloring Benchmarks
- CG, containing c++ codes for column generation BP and python codes for analyzing the results.- BP, contain the c++ code for branch-and-price and python codes for analyzing the results.

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

## Reproduce the experimental results in three steps
1. python3 01-train-and-optimize.py
2. python3 02-cg.py
3. python3 03-bp.py

For the second and third step, you can specificy the number of available CPUs (nCPUs).

## Results
There are two folders containing all the numerical results and figures in the main paper and the appendix:
- `results_cg' contains the results for column generation
- `results_bp' containing the results for branch-and-price

The contents of a result file is a latex table or data that can be ploted (using tikz) in latex directly.

## Assiciation between result files and Tables/Figures in the paper:
- 
