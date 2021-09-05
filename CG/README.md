# A Machine-Learning-Based Pricing Heuristic for Column Generation

argv[1] options:

- 0: train SVM
- 1: train Logistic Model           
    - argv[2]: b0
    - argv[3]: b1

- 4-11: test heuristic pricing methods 
    - argv[1]:
        - 6: test mlph
        - 7: test aco 
        - 8: test tsm 
        - 9: test lscc 
        - 10: test fastwclq 
        - 11: test gurobi  
    - argv[2] is the benchmark set: 0 for small-scale graphs and 1 for large-scale graphs
    - argv[3] is the test problem index, and -1 for test all problems
    - argv[4] is the column selection method
        - 0: add-partial
        - 1: add-all
        - 5: replace-existing
    - argv[5] is the seed for generating initial random columns
