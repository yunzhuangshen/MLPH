from utils import *

def find_insts_solved_by_all_methods(seeds, methods, method_rets, inst_names):

    insts = []
    for inst_id, inst_name in enumerate(inst_names): 
        tot_optimality = 0
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

                optimality = int(optimality)
                tot_optimality += optimality
        if tot_optimality == len(seeds) * len(methods):
            insts.append(inst_name)
    return insts


def make_table(seeds, methods, method_rets, inst_names, outdir_name):

    ret = [[ np.zeros(4) for j in range(len(methods))] for i in range(len(inst_names))]

    for inst_id, inst_name in enumerate(inst_names): 
        for method_id, (method, method_ret) in enumerate(zip(methods, method_rets)):
            for seed in seeds:
                subdir_path = join(method_ret, f'seed_{seed}')
                solving_fpath = join(subdir_path, f'{inst_name}_solving_stats.csv')
                cg_fpath = join(subdir_path, f'{inst_name}_cg_stats.csv')
                if not exists(solving_fpath):
                    print(f'solving stats file not exist: {solving_fpath}')
                    assert(False)
                with open (solving_fpath,'r') as solving_file:
                    tokens = solving_file.readlines()[1].strip().split(',')
                isOptimal, cg_iters, heur_successes = int(tokens[0]), int(tokens[7]), int(tokens[9])
                cg_iters = int (cg_iters)
                print(solving_fpath, cg_iters)
                with open (cg_fpath,'r') as cg_file:
                    first_iter_tokens = cg_file.readlines()[1].strip().split(',')
                nb_nrc_cols = float(first_iter_tokens[3])
                min_nrc = float(first_iter_tokens[4])
                mean_nrc=float(first_iter_tokens[5])
                median=float(first_iter_tokens[6])
                stdev_nrc=float(first_iter_tokens[7])

                tmp=np.array([cg_iters ,nb_nrc_cols, min_nrc, mean_nrc])
                ret[inst_id][method_id] += tmp

            ret[inst_id][method_id][0] /= len(seeds)
            ret[inst_id][method_id][1] /= len(seeds)
            ret[inst_id][method_id][2] /= len(seeds)
            ret[inst_id][method_id][3] /= len(seeds)


    stats_fname = f"{outdir_name}/heuristic_cg_relation.csv"
    stats_file = open(stats_fname, 'w')
    for inst_id, inst_name in enumerate(inst_names): 
        stats_file.write(inst_name)
        for val_id in range(3):
            for method_id, method in enumerate(methods):
                tmp = round(ret[inst_id][method_id][val_id],2)
                stats_file.write(f',{tmp}')
        stats_file.write('\n')

    stats_file.close()
    os.system(f'tably -n {stats_fname} > {outdir_name}/heuristic_cg_relation.tex')
    os.remove(stats_fname)

if __name__ == '__main__':
    seeds = [i for i in range(1, 25)]
    # methods = ['svm_cs_0','aco','gurobi','gurobi_heur','tsm','fastwclq','lscc']
    methods = ['svm_cs_0']
    method_rets = [f'results_small/results_{m}' for m in methods]
    inst_names, t_means = get_test_inst_names_lp_sorted(seeds, methods, method_rets, 'small')
    os.makedirs('data/small',exist_ok=True)
    inst_names = find_insts_solved_by_all_methods(seeds, methods, method_rets, inst_names)
    inst_names = ['4-Insertions_3']
    make_table(seeds, methods, method_rets, inst_names, 'data/small')

