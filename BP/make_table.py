from utils import *

# column names
OPTIMALITY = 0; ROOT_SOLVED=1; GAP = 2; TOTAL_TIME = 3
NNODES = 4; ROOT_TIME = 5; EXACT_TIME = 6; HEUR_TIME = 7; EXACT_AND_HEUR_TIME = 8

def gen_latex(records_all, columns, inst_names, methods, winning_method_name, prefix_path, dest_dir):

    gstats_dic = get_gstats(prefix_path)
    os.makedirs(dest_dir, exist_ok=True)
    stats_fname = f'{dest_dir}/win_{winning_method_name}.csv'
    stats_file = open(stats_fname, 'w')

    for inst_name in inst_names:
        record = records_all[inst_name]
        nb_node, density = gstats_dic[inst_name]
        stats_file.write(f',{inst_name},{nb_node},{density}')
        if OPTIMALITY in columns:
            for method_id, method in enumerate(methods):
                stats_file.write(f",{int(record[method_id][OPTIMALITY])}")
        if ROOT_SOLVED in columns:
            for method_id, method in enumerate(methods):
                stats_file.write(f",{int(record[method_id][ROOT_SOLVED])}")
        if GAP in columns:
            for method_id, method in enumerate(methods):
                if record[method_id][GAP] == -1: #placeholder for N/A
                    stats_file.write(f",N/A")
                else:
                    stats_file.write(f",{round(record[method_id][GAP], 1)}")
        if TOTAL_TIME in columns:
            for method_id, method in enumerate(methods):
                stats_file.write(f",{int(record[method_id][TOTAL_TIME])}")
        if ROOT_TIME in columns:
            for method_id, method in enumerate(methods):
                stats_file.write(f",{int(record[method_id][ROOT_TIME])}")
        if HEUR_TIME in columns:
            for method_id, method in enumerate(methods):
                stats_file.write(f",{int(record[method_id][HEUR_TIME])}")
        if EXACT_TIME in columns:
            for method_id, method in enumerate(methods):
                stats_file.write(f",{int(record[method_id][EXACT_TIME])}")
        if NNODES in columns:
            for method_id, method in enumerate(methods):
                stats_file.write(f",{int(record[method_id][NNODES])}")
        stats_file.write('\n')
    stats_file.close()
    os.system(f'tably -n {stats_fname} > {dest_dir}/win_{winning_method_name}.tex')
    os.remove(stats_fname)

def classify_insts_based_on_comparable_results(methods, method_rets, inst_names, seeds):

    ret = [[ 0 for j in range(len(methods))] for i in range(len(inst_names))]

    insts_solved_by_both = []
    insts_not_solved_by_both = []
    insts_the_rest = []

    for inst_id, inst_name in enumerate(inst_names): 
        for method_id, (method, method_ret) in enumerate(zip(methods, method_rets)):
            for seed in seeds:
                ret_path = f'{method_ret}/{seed}/{inst_name}.txt'
                with open(ret_path, 'r') as f:
                    optimality, primal_bound, dual_bound, gap, nnodes, tot_time, root_time, exact_time, heur_time, fixing_col_time, primal_heur_time = f.readlines()[0].strip().split(',')
                    
                # OPTIMALITY = 0, ROOT_SOLVED=1, GAP = 2, TOTAL_TIME = 3, NNODES = 4, ROOT_TIME = 5, 
                # ROOT_SOLVED = 6, EXACT_TIME = 7, HEUR_TIME = 8
                
                optimality = int(optimality)
                ret[inst_id][method_id] += optimality

        if ret[inst_id][0] + ret[inst_id][1] == 2*len(seeds):
            insts_solved_by_both.append(inst_name)
        elif ret[inst_id][0] + ret[inst_id][1] == 0:
            insts_not_solved_by_both.append(inst_name)
        else:
            insts_the_rest.append(inst_name)

    return insts_solved_by_both, insts_not_solved_by_both, insts_the_rest

def make_stats_table(methods, method_rets, inst_names, seeds, win_dict, prefix_path):

    records_all = {}
    ret = [[ np.zeros(8) for j in range(len(methods))] for i in range(len(inst_names))]

    for inst_id, inst_name in enumerate(inst_names): 
        for method_id, (method, method_ret) in enumerate(zip(methods, method_rets)):
            for seed in seeds:
                ret_path = f'{method_ret}/{seed}/{inst_name}.txt'
                with open(ret_path, 'r') as f:
                    optimality, primal_bound, dual_bound, gap, nnodes, tot_time, root_time, exact_time, heur_time, fixing_col_time, primal_heur_time = f.readlines()[0].strip().split(',')
                    
                # OPTIMALITY = 0, ROOT_SOLVED=1, GAP = 2, TOTAL_TIME = 3, NNODES = 4, ROOT_TIME = 5, 
                # ROOT_SOLVED = 6, EXACT_TIME = 7, HEUR_TIME = 8
                
                optimality = int(optimality)
                ret[inst_id][method_id][OPTIMALITY] += optimality

                tot_time = float(tot_time)
                dual_bound = ceil(float(dual_bound))
                primal_bound = float(primal_bound)
                nnodes = float(nnodes)
                root_time = float(root_time)
                exact_time = float(exact_time)
                heur_time = float(heur_time)
                solve_root = dual_bound != -1e20
                if solve_root:
                    ret[inst_id][method_id][ROOT_SOLVED] += 1
                    gap = (primal_bound - dual_bound)/primal_bound
                    ret[inst_id][method_id][GAP] += gap*100

                ret[inst_id][method_id][TOTAL_TIME] += tot_time
                ret[inst_id][method_id][NNODES] += nnodes
                ret[inst_id][method_id][ROOT_TIME] += root_time
                ret[inst_id][method_id][EXACT_TIME] += exact_time
                ret[inst_id][method_id][HEUR_TIME] += heur_time
            
            # OPTIMALITY = 0, ROOT_SOLVED=1, GAP = 2, TOTAL_TIME = 3, 
            # NNODES = 4, ROOT_TIME = 5, EXACT_TIME = 6, HEUR_TIME = 7
            
            if ret[inst_id][method_id][ROOT_SOLVED] > 0:
                ret[inst_id][method_id][GAP]/=ret[inst_id][method_id][ROOT_SOLVED]
            else:
                ret[inst_id][method_id][GAP]=-1

            ret[inst_id][method_id][3:]/=len(seeds)
        records_all[inst_name] = ret[inst_id]


    dest_dir_prefix = f'{prefix_path}/../results_bp/table_{methods[0]}-{methods[1]}'
    dest_dir = f'{dest_dir_prefix}/full/'
    os.makedirs(dest_dir, exist_ok=True)
    method1_win_insts = [item[0] for item in win_dict[methods[0]]]
    method2_win_insts = [item[0] for item in win_dict[methods[1]]]
    none_win_insts = win_dict['none']
    
    # generate table containing full statistics
    # columns = [OPTIMALITY, ROOT_SOLVED, GAP, TOTAL_TIME,HEUR_TIME,EXACT_TIME,NNODES]
    columns = [OPTIMALITY, ROOT_SOLVED, GAP, TOTAL_TIME]
    gen_latex(records_all, columns, method1_win_insts, methods, methods[0], prefix_path, dest_dir)
    gen_latex(records_all, columns, method2_win_insts, methods, methods[1], prefix_path, dest_dir)
    gen_latex(records_all, columns, none_win_insts, methods, 'none', prefix_path, dest_dir)

    insts_solved_by_both, insts_not_solved_by_both, insts_the_rest = \
        classify_insts_based_on_comparable_results(methods, method_rets, inst_names, seeds)
    
    # generate table for computational time for instances solved by both methods
    dest_dir = f'{dest_dir_prefix}/time_for_all_solved/'
    os.makedirs(dest_dir, exist_ok=True)
    gen_latex(records_all, [TOTAL_TIME], [inst for inst in method1_win_insts if inst in insts_solved_by_both], methods, methods[0], prefix_path, dest_dir)
    gen_latex(records_all, [TOTAL_TIME], [inst for inst in method2_win_insts if inst in insts_solved_by_both], methods, methods[1], prefix_path, dest_dir)
    gen_latex(records_all, [TOTAL_TIME], [inst for inst in none_win_insts if inst in insts_solved_by_both], methods, 'none', prefix_path, dest_dir)

    # generate table for rootsolve and gap for instances not solved by either of the two methods
    dest_dir = f'{dest_dir_prefix}/gap_for_all_not_solved/'
    os.makedirs(dest_dir, exist_ok=True)
    gen_latex(records_all, [ROOT_SOLVED, GAP], [inst for inst in method1_win_insts if inst in insts_not_solved_by_both], methods, methods[0], prefix_path, dest_dir)
    gen_latex(records_all, [ROOT_SOLVED, GAP], [inst for inst in method2_win_insts if inst in insts_not_solved_by_both], methods, methods[1], prefix_path, dest_dir)
    gen_latex(records_all, [ROOT_SOLVED, GAP], [inst for inst in none_win_insts if inst in insts_not_solved_by_both], methods, 'none', prefix_path, dest_dir)

    
    # generate table for number of solves for the rest of the instances
    dest_dir = f'{dest_dir_prefix}/number_solve_for_not_all_solved'
    os.makedirs(dest_dir, exist_ok=True)
    gen_latex(records_all, [OPTIMALITY], [inst for inst in method1_win_insts if inst in insts_the_rest], methods, methods[0], prefix_path, dest_dir)
    gen_latex(records_all, [OPTIMALITY], [inst for inst in method2_win_insts if inst in insts_the_rest], methods, methods[1], prefix_path, dest_dir)
    gen_latex(records_all, [OPTIMALITY], [inst for inst in none_win_insts if inst in insts_the_rest], methods, 'none', prefix_path, dest_dir)


def make_tables(method1, method2, prefix_path):
    seeds = [i for i in range(1, 25)]
    methods = [method1, method2]; method_rets = [f'{prefix_path}/results/{m}' for m in methods]
    # unknown bug to 'mug' graphs, resulting in solving to optimality at cutoff time
    inst_names = get_test_inst_names(methods[0], seeds[0], prefix_path)
    inst_names = [inst_name for  inst_name in inst_names if 'mug' not in inst_name]

    win_dict = compare(methods, method_rets, inst_names, seeds, prefix_path)
    # for k,v in win_dict.items():
    #     print(k, len(v))
    make_stats_table(methods, method_rets, inst_names, seeds, win_dict, prefix_path)
