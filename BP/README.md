# A Machine-Learning-Based Pricing Heuristic for Column Generation

command to run: ./BP [graph name] [method] [seed] 

Method:
- 0 BP-DEF
- 1 BP-None
- 2 BP-MLPH (several parameters to specifiy following the parameter [seed])
    - argv[4]: sample_size
    - argv[5]: root_node_column_selection_limit
    - argv[6]: child_node_column_selection_limit
- 3 BP-MLPH-force-exact

