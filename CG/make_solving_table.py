from utils import *

def make_solving_table(seeds, methods, method_rets, inst_names, outdir_name):

    ret = [[ np.zeros(3) for j in range(len(methods))] for i in range(len(inst_names))]
    best_lps = np.zeros(len(inst_names))
    best_lps_optimal = np.array([False]*len(inst_names))
    for inst_id, inst_name in enumerate(inst_names): 
        best_lp_obj = 1e8; best_lp_optimal=False
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

                current_optimal = optimality=='1'
                current_lp_obj = float(lp_obj)

                # sanity check, optimal lp obj should be the same
                if current_optimal and best_lp_optimal:
                    # if '4-FullIns_5' == inst_name:
                    #     print(method, seed)

                    if best_lp_obj != current_lp_obj:
                        print(solving_stats_path)
                        print(inst_name, best_lp_obj, current_lp_obj)
                    assert(best_lp_obj == current_lp_obj)
                if current_lp_obj < best_lp_obj:
                    best_lp_obj = current_lp_obj
                best_lp_optimal = best_lp_optimal or current_optimal

                # log statistics
                opt = float(tot_time) < 8100 and int(optimality) == 1
                ret[inst_id][method_id][0]+=int(opt)
                ret[inst_id][method_id][1]+=float(lp_obj)
                ret[inst_id][method_id][2]+=float(tot_time)
            ret[inst_id][method_id][1:] /= len(seeds)
        best_lps[inst_id] = best_lp_obj
        best_lps_optimal[inst_id] = best_lp_optimal

    stats_fname = f"{outdir_name}/solving_stats.csv"
    stats_file = open(stats_fname, 'w')
    for inst_id, inst_name in enumerate(inst_names): 
        
        stats_file.write(inst_name)
        for val_id in range(3):
            for method_id, method in enumerate(methods):
                tmp = 'N/A'
                if val_id==0:
                    tmp = int(ret[inst_id][method_id][val_id])
                elif val_id==1:
                    tmp = round(ret[inst_id][method_id][val_id],3)
                else:
                    tmp = round(ret[inst_id][method_id][val_id],1)

                stats_file.write(f',{tmp}')
        stats_file.write('\n')

    stats_file.close()
    with open(stats_fname, 'r') as f:
        M = genfromtxt(stats_fname, delimiter=',')[:, 1:]
    nmethod = len(methods)

    # easy_insts = np.array(inst_names)[np.sum(M[:, :nmethod], axis=1) > 0].tolist()
    # hard_insts = np.array(inst_names)[np.sum(M[:, :nmethod], axis=1) == 0].tolist()
    
    # print(f'easy_insts: {len(easy_insts)}')
    # print(sorted(easy_insts))
    # print(f'hard_insts: {len(hard_insts)}')
    # print(sorted(hard_insts))

    tmp=[]
    for col_idx in range(nmethod * 3):
        if col_idx < nmethod:
            tmp.append(str(int(np.sum(M[:,col_idx]).item())))
        elif col_idx < nmethod*2:
            tmp.append(str(round(gmean(M[:,col_idx]).item(),3)))
        else:
            tmp.append(str(round(gmean(M[:,col_idx]).item(),1)))

    
    with open(stats_fname, 'a') as f:
        f.write('statistic,' + ','.join(tmp) + '\n')
    os.system(f'tably -n {stats_fname} > {outdir_name}/table_solving_stats.tex')
    os.remove(stats_fname)