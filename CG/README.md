# A Machine-Learning-Based Pricing Heuristic for Column Generation

argv[1] options:

- 0: train SVM
- 1: train Logistic Model           
    - argv[2]: b0
    - argv[3]: b1

- 4-11: test heuristic pricing methods 
    - argv[1]:
        - 6: test MLPH
        - 7: test ACO 
        - 8: test TSM 
        - 9: test LSCC 
        - 10: test Fastwclq 
        - 11: test Gurobi  
    - argv[2] is the benchmark set: 0 for small-scale graphs and 1 for large-scale graphs
    - argv[3] is the test problem index. -1 for test all problems
    - argv[4] is the column selection method
        - 0: add-partial
        - 1: add-all
        - 5: replace-existing
    - argv[5] is the seed for generating initial random columns
