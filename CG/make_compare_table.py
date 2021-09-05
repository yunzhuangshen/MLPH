from utils import *

def find_nth_place(inst_name, seeds, methods, method_rets, place):
    ret = [ np.zeros(2) for j in range(len(methods))]

    for method_id, (method, method_ret) in enumerate(zip(methods, method_rets)):
        for seed in seeds:
            solving_stats_path = join(method_ret, f'seed_{seed}', f'{inst_name}_solving_stats.csv')
            if not exists(solving_stats_path):
                print(f'solving stats file not exist: {solving_stats_path}')
                assert(False)
            with open(solving_stats_path, 'r') as f:  
                lines =  f.readlines()
            if len(lines) <= 1:
                print(f'INCOMPLETE: {solving_stats_path}')
                assert(False)
            tokens = lines[1].strip().split(',')
            optimality, lp_obj, tot_time, cpu_time, master_duration,\
            heur_pricing_duration, exact_pricing_duration,\
            num_cg_iter, num_added_cols, num_heur_success = tokens[:10]    

            # log statistics
            ret[method_id][0]+=float(lp_obj)
            ret[method_id][1]+=float(tot_time)
        ret[method_id] /= len(seeds)
    
    masked_method_ids=[]
    while place != 0:
        best_lp_obj = 1e8; wining_method_id = -1
        for method_id, method in enumerate(methods):
            if method_id not in masked_method_ids:
                if best_lp_obj>ret[method_id][0]:
                    best_lp_obj=ret[method_id][0]

        win_methods_ids = []
        for method_id, method in enumerate(methods):
            if method_id not in masked_method_ids:
                if best_lp_obj==ret[method_id][0]:
                    win_methods_ids.append(method_id)
        
        if len(win_methods_ids) == 1:
            wining_method_id = win_methods_ids[0]
        else:
            shorest_time = 1e8
            for method_id in win_methods_ids:
                if shorest_time>ret[method_id][1]:
                    shorest_time=ret[method_id][1]
                    wining_method_id=method_id
        
        place-=1
        masked_method_ids.append(wining_method_id)
    return wining_method_id, methods[wining_method_id]

def make_compare_table(prefix='./', place=1):
    seeds = [i for i in range(1, 25)]
    methods = ['mlph_cs0','aco','gurobi','gurobi_heur','tsm','fastwclq','lscc']
    method_rets = [f'{prefix}/results_small/{m}' for m in methods]
    g_names, _ = get_test_inst_names_lp_sorted(seeds, methods, method_rets, 'small', prefix=prefix)
    g_paths = [f'{prefix}/../GCB/{g_name}.col' for g_name in g_names]
    ret=[]
    count = [0 for _ in range(len(methods))]
    for g_name, g_path in zip(g_names, g_paths):
        with open(g_path) as f:
            lines = f.readlines()
        for line in lines:
            if line[0] == 'p':
                break
        _, _, nb_node, nb_edge = line.split()
        nb_node, nb_edge = int(nb_node), int(nb_edge)
        density = 2 * nb_edge / (nb_node * (nb_node-1))
        if density>1:
            density = density / 2
        wining_method_id, wining_method = find_nth_place(g_name, seeds, methods, method_rets, place)
        count[wining_method_id]+=1
        ret.append(f'{g_name},{nb_node},{round(density, 3)},{wining_method}')

    with open(f'{prefix}/../results_cg/compare_number.txt', 'w') as f:
        for method, ctr in zip(methods, count):
            f.write(f'{method},{ctr}\n')

    with open(f'{prefix}/../results_cg/compare_figure.txt', 'w') as f:
        f.write('x,y,label\n')
        for line in ret:
            _, nb_node, density, wining_method = line.split(',')
            nb_node = float(nb_node)
            density = float(density)
            f.write(f'{nb_node},{density},{wining_method}\n')