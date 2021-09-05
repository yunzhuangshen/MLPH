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
        # optimal_count = 0
        # for method_id, method in enumerate(methods):
        #     optimal_count += ret[inst_id][method_id][0]

        # if optimal_count == 0:
        #     continue
        
        stats_file.write(inst_name)

        # avg cg iter success, avg cg iter unsuccess, success_rate.
        for val_id in range(1,4):
            for method_id, method in enumerate(methods):
                tmp = 'N/A'
                if (val_id==1 or val_id==2):
                    if not np.isnan(ret[inst_id][method_id][val_id]):
                        tmp = round(ret[inst_id][method_id][val_id],1)
                elif val_id==3:
                    tmp = round(ret[inst_id][method_id][val_id],2)
                else:
                    tmp = round(ret[inst_id][method_id][val_id],1)
                
                stats_file.write(f',{tmp}')
        stats_file.write('\n')

    stats_file.close()
    os.system(f'tably -n {stats_fname} > {outdir_name}/table_cg.tex')
    os.remove(stats_fname)