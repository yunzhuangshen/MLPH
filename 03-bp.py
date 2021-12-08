import sys, os, time, subprocess
sys.path.append('BP/')
from math import floor
from BP.analyze import analyze,analyze_MLPH_variants

BP_DEF=0; BP_None=1; BP_MLPH=2; BP_MLPH_FORCE_EXACT=3
seeds = [i for i in range(1, 25)]
insts = ['1-FullIns_4', '1-FullIns_5', '1-Insertions_4', '1-Insertions_5', '2-FullIns_4', '2-FullIns_5', '2-Insertions_3', '2-Insertions_4', '2-Insertions_5', '3-FullIns_4', '3-FullIns_5', '3-Insertions_3', '3-Insertions_4', '3-Insertions_5', '4-FullIns_4', '4-FullIns_5', '4-Insertions_3', '4-Insertions_4', '5-FullIns_4', 'DSJC1000.9', 'DSJC125.1', 'DSJC125.5', 'DSJC125.9', 'DSJC250.5', 'DSJC250.9', 'DSJC500.9', 'DSJR500.1c', 'DSJR500.5', 'ash331GPIA', 'ash608GPIA', 'flat300_20_0', 'flat300_26_0', 'flat300_28_0', 'latin_square_10', 'le450_15a', 'le450_15b', 'le450_25a', 'le450_25b', 'le450_25c', 'le450_25d', 'le450_5c', 'le450_5d', 'myciel5', 'myciel6', 'myciel7', 'qg.order100', 'qg.order30', 'qg.order40', 'qg.order60', 'queen10_10', 'queen11_11', 'queen12_12', 'queen13_13', 'queen14_14', 'queen15_15', 'queen16_16', 'queen9_9', 'r1000.1c', 'r1000.5', 'r125.5', 'r250.5', 'school1', 'school1_nsh', 'wap05a', 'wap06a', 'will199GPIA'] 
sample_factor=10.; root_cs_limit_factor=1.; child_cs_limit_factor=0.1

# insts = ['le450_25a'] 
# seeds = [1,2]

def test():
    test_methods = [BP_DEF, BP_MLPH]

    triples = [(test_method, inst, seed) for test_method in test_methods for inst in insts for seed in seeds]
    nb_runs = len(triples)
    print(f'total number of runs: {nb_runs}')

    triples = (trip for trip in triples)
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
        while len(proc_pool) < nCPUs and hasNext:
            item = next(triples,None)
            hasNext = item is not None
            if hasNext:
                (method, inst, seed) = item
                proc = subprocess.Popen(f'./BP {inst} {method} {seed} {sample_factor} {root_cs_limit_factor} {child_cs_limit_factor} > {method}_{inst}_{seed}_{sample_factor}_{root_cs_limit_factor}_{child_cs_limit_factor}.log', cwd=f'./BP/build', shell=True)
                proc_pool.append(proc)
            
        # report progress
        curr_time = time.time()
        cur_duration = floor((curr_time - start_time)/3600)
        if (cur_duration > duration):
            duration = cur_duration
            print(f"time used: {duration} hours\nnumber of finished runs: {nb_finished}\nnumber of remaining runs: {nb_runs - nb_finished}\n")
        time.sleep(5)

def test_variants():

    # BP without any heuristic pricing method
    quadruples = [(BP_None, inst, seed, None) for inst in insts for seed in seeds]
    # BP-MLPH with forcing exact method for overcoming tailoff effect
    quadruples.extend([(BP_MLPH_FORCE_EXACT, inst, seed, None) for inst in insts for seed in seeds])
    # BP-MLPH with other parameters
    quadruples.extend([(BP_MLPH, inst, seed, (10,1,1)) for inst in insts for seed in seeds])
    quadruples.extend([(BP_MLPH, inst, seed, (1,1,0.1)) for inst in insts for seed in seeds])
    quadruples.extend([(BP_MLPH, inst, seed, (0.1,1,0.1)) for inst in insts for seed in seeds])

    
    nb_runs = len(quadruples)
    print(f'total number of runs: {nb_runs}')

    quadruples = (quad for quad in quadruples)
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
        while len(proc_pool) < nCPUs and hasNext:
            item = next(quadruples,None)
            hasNext = item is not None
            if hasNext:
                (method, inst, seed, params) = item
                if method == BP_MLPH:
                    proc = subprocess.Popen(f'./BP {inst} {method} {seed} {params[0]} {params[1]} {params[2]} > {method}_{inst}_{seed}_{params[0]}_{params[1]}_{params[2]}.log', cwd=f'./BP/build', shell=True)
                else:
                    proc = subprocess.Popen(f'./BP {inst} {method} {seed} > {method}_{inst}_{seed}.log', cwd=f'./BP/build', shell=True)
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
    nCPUs = 4; 
    os.system(f'mkdir BP/build; cd BP/build && cmake ../ && make')
    test(); test_variants(); 
    analyze(prefix='./BP/')
    analyze_MLPH_variants(prefix='./BP/')