import sys, os, time, subprocess
sys.path.append('CG/')
from math import floor
from CG.analyze import analyze

SMALL_GRAPHS=0; LARGE_GRAPHS=1
MLPH=6; ACO=7; TSM=8; LSCC=9; FASTWCLQ=10; GUROBI=11; GUROBI_HEUR=12
ADD_PARTIAL=0; ADD_ALL=1; REPLACE_EXISTING=5
seeds = [i for i in range(1, 25)]

def test_pricing_methods(graph_type):
    global nCPUs
    if graph_type == SMALL_GRAPHS:
        test_methods = [MLPH, ACO, TSM, LSCC, FASTWCLQ, GUROBI, GUROBI_HEUR]
        nProcesses = nCPUs
    else:
        test_methods = [MLPH, ACO, FASTWCLQ, GUROBI, GUROBI_HEUR]
        nProcesses = max(1, nCPUs//4) # large graphs require 4 CPUs each run

    pairs = [(test_method, seed) for test_method in test_methods for seed in seeds]
    nb_runs = len(pairs)
    print(f'total number of runs for {"small" if graph_type == SMALL_GRAPHS else "large"} graphs: {nb_runs}')
    print(f'each run solves the GCP-LPs on the set of {81 if graph_type == SMALL_GRAPHS else 8} graphs')

    pairs = (p for p in pairs)
    proc_pool = []
    start_time = time.time()
    nb_finished = 0    
    duration = 0
    hasNext = True

    while hasNext or len(proc_pool) > 0:
        
        # remove fininshed processes
        for idx in reversed(range(len(proc_pool))):
            if proc_pool[idx].poll() is not None:
                del proc_pool[idx]
                nb_finished+=1
        
        # execute new runs
        while len(proc_pool) < nProcesses and hasNext:
            item = next(pairs,None)
            hasNext = item is not None
            if hasNext:
                method, seed = item
                proc = subprocess.Popen(f'./CG {method} {graph_type} -1 {ADD_PARTIAL} {seed} > g{graph_type}_m{method}_cs{ADD_PARTIAL}_s{seed}.log', cwd=f'./CG/build', shell=True)
                proc_pool.append(proc)
            
        # report progress
        curr_time = time.time()
        cur_duration = floor((curr_time - start_time)/3600)
        if (cur_duration > duration):
            duration = cur_duration
            print(f"time used: {duration} hours\nnumber of finished runs: {nb_finished}\nnumber of remaining runs: {nb_runs - nb_finished}\n")
        time.sleep(5)

        

def test_column_selection():
    global nCPUs

    nProcesses = max(1, nCPUs//4) # large graphs require 4 CPUs each run
    cs_methods = [ADD_ALL, REPLACE_EXISTING]
    pairs = [(test_method, seed) for test_method in cs_methods for seed in seeds]
    nb_runs = len(pairs)
    print(f'total number of runs for column selection methods: {nb_runs}')

    pairs = (p for p in pairs)
    proc_pool = []
    start_time = time.time()
    nb_finished = 0    
    duration = 0
    hasNext = True

    while hasNext or len(proc_pool) > 0:
        
        # remove fininshed processes
        for idx in reversed(range(len(proc_pool))):
            if proc_pool[idx].poll() is not None:
                del proc_pool[idx]
                nb_finished+=1

        # execute new runs
        while len(proc_pool) < nProcesses and hasNext:
            item = next(pairs,None)
            hasNext = item is not None
            if hasNext:
                method, seed = item
                proc = subprocess.Popen(f'./CG {MLPH} {LARGE_GRAPHS} -1 {method} {seed} > g{LARGE_GRAPHS}_m{MLPH}_cs{method}_s{seed}.log', cwd=f'./CG/build', shell=True)
                proc_pool.append(proc)

        # report progress
        curr_time = time.time()
        cur_duration = floor((curr_time - start_time)/3600)
        if (cur_duration > duration):
            duration = cur_duration
            print(f"time used: {duration} hours\nnumber of finished runs: {nb_finished}\nnumber of remaining runs: {nb_runs - nb_finished}\n")
        time.sleep(5)

if __name__ =='__main__':
    
    global nCPUs
    nCPUs = 4; # should be multiples of 4 for testing on large graphs

    os.system(f'mkdir CG/build; cd CG/build && cmake ../ && make')
    test_pricing_methods(SMALL_GRAPHS)
    test_pricing_methods(LARGE_GRAPHS)
    test_column_selection()
    analyze(prefix='./CG/')

