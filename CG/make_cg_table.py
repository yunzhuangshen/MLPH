from utils import *

def make_cg_table(seeds, methods, method_rets, inst_names, outdir_name):

    ret = [[ np.zeros(6) for j in range(len(methods))] for i in range(len(inst_names))]

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
                isOptimal, master_duration, pricing_duration, exact_pricing_duration, cg_iters, heur_successes = \
                        int(tokens[0]), float(tokens[4]), float(tokens[5]), float(tokens[6]), int(tokens[7]), int(tokens[9])
                
                divider = cg_iters-1 if isOptimal==1 else cg_iters
                if divider==0:
                    divider=1
                tmp=np.array([isOptimal, cg_iters if isOptimal==1 else 0, cg_iters if isOptimal==0 else 1,\
                    heur_successes/divider, master_duration, pricing_duration+exact_pricing_duration])

                ret[inst_id][method_id] += tmp
        
            nsuccess = ret[inst_id][method_id][0]
            nfail = len(seeds) - nsuccess
            if nsuccess == 0:
                ret[inst_id][method_id][1] = np.nan
            else:
                ret[inst_id][method_id][1] /= nsuccess 
            if nfail==0:
                ret[inst_id][method_id][2] = np.nan
            else:
                ret[inst_id][method_id][2] /= nfail

            ret[inst_id][method_id][3] /= len(seeds)
            ret[inst_id][method_id][4] /= len(seeds)
            ret[inst_id][method_id][5] /= len(seeds)


    stats_fname = f"{outdir_name}/cg_stats.csv"
    stats_file = open(stats_fname, 'w')
    for inst_id, inst_name in enumerate(inst_names): 
        optimal_count = 0
        for method_id, method in enumerate(methods):
            optimal_count += ret[inst_id][method_id][0]

        if optimal_count == 0:
            continue
        
                # avg cg iter success, avg cg iter unsuccess, success_rate.

        stats_file.write(f'{inst_name}')
        for method_id, method in enumerate(methods):
            tmp1 = int(ret[inst_id][method_id][0])
            tmp2 = 'N/A' if np.isnan(ret[inst_id][method_id][1]) else round(ret[inst_id][method_id][1],1)
            stats_file.write(f',{tmp2} ({tmp1})')
        for method_id, method in enumerate(methods):
            tmp1 = len(seeds) - int(ret[inst_id][method_id][0])
            tmp2 = 'N/A' if np.isnan(ret[inst_id][method_id][2]) else round(ret[inst_id][method_id][2],1)
            stats_file.write(f',{tmp2} ({tmp1})')
        stats_file.write(f'\n')

    stats_file.close()
    os.system(f'tably -n {stats_fname} > {outdir_name}/table_cg.tex')
    os.remove(stats_fname)

if __name__ == '__main__':
    prefix='./'
    seeds = [i for i in range(1, 25)]
    methods = ['mlph_cs0','aco','gurobi','gurobi_heur','tsm','fastwclq','lscc']
    method_rets = [f'{prefix}/results_small/{m}' for m in methods]
    inst_names, t_means = get_test_inst_names_lp_sorted(seeds, methods, method_rets, 'small', prefix=prefix)
    dest_dir = f'{prefix}/../results_cg/small/'
    os.makedirs(dest_dir,exist_ok=True)
    make_cg_table(seeds, methods, method_rets, inst_names, dest_dir)