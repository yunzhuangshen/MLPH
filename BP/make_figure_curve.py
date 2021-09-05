from utils import *

def solving_time_curve(method1, method2, prefix_path='./', interval=10):
    
    seeds = [i for i in range(1,25)]
    methods = [method1, method2]; method_rets = [f'{prefix_path}/results/{m}' for m in methods]
    # unknown bug to 'mug' graphs, resulting in solving to optimality at cutoff time
    inst_names = get_test_inst_names(methods[0], seeds[0], prefix_path)
    inst_names = [inst_name for  inst_name in inst_names if 'mug' not in inst_name]
    
    
    seeds = [i for i in range(1,25)]
    plots = {}
    tot_time = 8000
    nslots = tot_time//interval
    for method in methods:
        plots[method] = [0 for i in range(nslots)]
            
    for method, method_ret_fname in zip(methods, method_rets):
        for inst_name in inst_names:
            for seed in seeds:
                subdir_path = join(method_ret_fname, f'{seed}')
                solving_stats_path = join(subdir_path, f'{inst_name}.txt')

                if not exists(solving_stats_path):
                    # 
                    # print(f'solving stats file not exist: {solving_stats_path}')
                    # assert(False)
                    pass
                else:
                    with open(solving_stats_path, 'r') as f:  
                        lines =  f.readlines()
                    if len(lines) < 1:
                        
                        print(f'ERROR at solving stats file: {solving_stats_path}')
                        assert(False)
                    
                    tokens = lines[0].strip().split(',')
                    optimal = tokens[0]=='1'
                    if optimal:
                        tokens = lines[0].strip().split(',')
                        idx = floor(float(tokens[5]) // interval)
                        for i in range(idx, nslots):
                            plots[method][i]+=1
    
    dest_dir = f'{prefix_path}/../results_bp/solving-curve'
    os.makedirs(dest_dir, exist_ok=True)
    for method, curve in plots.items():
        xs = [interval*i for i in range(1, nslots+1)]
        ys = curve
        assert(len(xs) == len(ys))
        with open(f'{dest_dir}/{method.replace("_","-")}.txt', 'w') as f:
            f.write('x,y\n')
            for i in range(len(xs)):
                f.write(f'{xs[i]/1000},{ys[i]}\n')

def solving_gap_curve(method1, method2, prefix_path, interval=1):
    seeds = [i for i in range(1,25)]
    methods = [method1, method2]; method_rets = [f'{prefix_path}/results/{m}' for m in methods]
    # unknown bug to 'mug' graphs, resulting in solving to optimality at cutoff time
    inst_names = get_test_inst_names(methods[0], seeds[0], prefix_path)
    inst_names = [inst_name for  inst_name in inst_names if 'mug' not in inst_name]
    
    max_gap = -1e8
    for method, method_ret_fname in zip(methods, method_rets):
        unsolved_seeds = {}
        bug_insts = {}; aborted_runs = {}
        for inst_name in inst_names:
            for seed in seeds:
                subdir_path = join(method_ret_fname, f'{seed}')
                solving_stats_path = join(subdir_path, f'{inst_name}.txt')

                if not exists(solving_stats_path):
                    
                    print(f'solving stats file not exist: {solving_stats_path}')
                    # assert(False)
                    if inst_name not in aborted_runs:
                        aborted_runs[inst_name] = [seed]
                    else:
                        aborted_runs[inst_name].append(seed)
                else:
                    with open(solving_stats_path, 'r') as f:  
                        lines =  f.readlines()
                    if len(lines) < 1:
                        
                        print(f'ERROR at solving stats file: {solving_stats_path}')
                        assert(False)
                    optimality, primal_bound, dual_bound, gap, nnodes, tot_time, root_time, exact_time, heur_time, fixing_col_time, primal_heur_time = lines[0].strip().split(',')
                    dual_bound = ceil(float(dual_bound))
                    primal_bound = float(primal_bound)
                    if optimality == '1' and float(tot_time) >= 8000:
                        if inst_name in bug_insts:
                            bug_insts[inst_name]+=1
                        else:
                            bug_insts[inst_name] = 1
                    solve_root = dual_bound != -1e20
                    if solve_root:
                        gap = (primal_bound - dual_bound)/primal_bound * 100
                        if max_gap < gap:
                            max_gap = gap
                    else:
                        if inst_name not in unsolved_seeds:
                            unsolved_seeds[inst_name] = []
                        unsolved_seeds[inst_name].append(seed)
        # print('method: ', method)
        # for k,v in bug_insts.items():
        #     print(k,v)
        # for k,v in unsolved_seeds.items():
        #     print(k, len(v))
        # print()
    # print('max_gap: ', max_gap)

    gap_dict = {}; unsolved_root_dict = {}
    for method in methods:
        gap_dict[method] = [0 for i in range(ceil(max_gap) + 1)]
        unsolved_root_dict[method] = 0

    for method, method_ret_fname in zip(methods, method_rets):
        for inst_name in inst_names:
            for seed in seeds:
                subdir_path = join(method_ret_fname, f'{seed}')
                solving_stats_path = join(subdir_path, f'{inst_name}.txt')

                if not exists(solving_stats_path):
                    # 
                    # print(f'solving stats file not exist: {solving_stats_path}')
                    # assert(False)
                    pass
                else:
                    with open(solving_stats_path, 'r') as f:  
                        lines =  f.readlines()
                    if len(lines) < 1:
                        
                        print(f'ERROR at solving stats file: {solving_stats_path}')
                        assert(False)
                    optimality, primal_bound, dual_bound, gap, nnodes, tot_time, root_time, exact_time, heur_time, fixing_col_time, primal_heur_time = lines[0].strip().split(',')
                    dual_bound = ceil(float(dual_bound))
                    primal_bound = float(primal_bound)

                    solve_root = dual_bound != -1e20
                    if solve_root:
                        gap = (primal_bound - dual_bound)/primal_bound * 100
                        gap_dict[method][ceil(gap)] += 1
                    else:
                        unsolved_root_dict[method]+=1

    for method in methods:
        gap_dict[method] = np.cumsum(gap_dict[method]).tolist()    

    dest_dir = f'{prefix_path}/../results_bp/gap_curve_{method1}-{method2}/'
    os.makedirs(dest_dir,exist_ok=True)
    for method, curve in gap_dict.items():
        ys = [int(item) for item in curve]
        xs = [interval*i for i in range(len(ys))]

        with open(f'{dest_dir}/{method.replace("_","-")}.txt', 'w') as f:
            f.write('x,y\n')
            for i in range(len(xs)):
                f.write(f'{xs[i]},{round(ys[i],5)}\n')

            f.write(f'\n nb unsolved root: {unsolved_root_dict[method]}')
    
