from utils import *

def make_rc_table(seeds, methods, method_rets, inst_names, outdir_name):

    ret = [[ np.zeros(3) for j in range(len(methods))] for i in range(len(inst_names))]
    geo_mean_nb_nrc_cols = [[] for method in methods]
    geo_mean_min_nrc = [[] for method in methods]
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
                
                with open (cg_fpath,'r') as cg_file:
                    first_iter_tokens = cg_file.readlines()[1].strip().split(',')
                nb_nrc_cols = float(first_iter_tokens[3])
                min_nrc = float(first_iter_tokens[4])
                mean_nrc=float(first_iter_tokens[5])
                median=float(first_iter_tokens[6])
                stdev_nrc=float(first_iter_tokens[7])

                tmp=np.array([nb_nrc_cols, min_nrc, mean_nrc])
                ret[inst_id][method_id] += tmp
                geo_mean_nb_nrc_cols[method_id].append(nb_nrc_cols)
                geo_mean_min_nrc[method_id].append(-min_nrc)
            ret[inst_id][method_id][0] /= len(seeds)
            ret[inst_id][method_id][1] /= len(seeds)
            ret[inst_id][method_id][2] /= len(seeds)

    stats_fname = f"{outdir_name}/rc_mean_stats.csv"
    stats_file = open(stats_fname, 'w')
    stats_file.write("method,num_nrc_col,min_nrc\n")
    for method_id, method in enumerate(methods):
        stats_file.write(f"{method},{round(gmean(geo_mean_nb_nrc_cols[method_id]),1)},{round(-gmean(geo_mean_min_nrc[method_id]),2)}\n")
    stats_file.close()

    stats_fname = f"{outdir_name}/rc_stats.csv"
    stats_file = open(stats_fname, 'w')
    for inst_id, inst_name in enumerate(inst_names): 
        # optimal_count = 0
        # for method_id, method in enumerate(methods):
        #     optimal_count += ret[inst_id][method_id][0]

        # if optimal_count == 0:
        #     continue
        
        stats_file.write(inst_name)
        for val_id in range(3):
            for method_id, method in enumerate(methods):
                if val_id==0:
                    tmp = round(ret[inst_id][method_id][val_id],1)
                else:
                    tmp = round(ret[inst_id][method_id][val_id],2)
                
                stats_file.write(f',{tmp}')
        stats_file.write('\n')

    stats_file.close()
    os.system(f'tably -n {stats_fname} > {outdir_name}/table_rc.tex')
    os.remove(stats_fname)

if __name__ == '__main__':
    prefix='./'
    seeds = [i for i in range(1, 25)]
    methods = ['mlph_cs0','aco','gurobi','gurobi_heur','tsm','fastwclq','lscc']
    method_rets = [f'{prefix}/results_small/{m}' for m in methods]
    inst_names, t_means = get_test_inst_names_lp_sorted(seeds, methods, method_rets, 'small', prefix=prefix)
    dest_dir = f'{prefix}/../results_cg/small/'
    os.makedirs(dest_dir,exist_ok=True)
    make_rc_table(seeds, methods, method_rets, inst_names, dest_dir)