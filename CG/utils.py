import sys, os, glob
import numpy as np
from os.path import join, isdir, isfile, exists
from os import listdir
from math import floor
from numpy import genfromtxt

def gmean(arr, mask=None):
    arr = np.array(arr, dtype=float)
    shift=1
    log_a = np.log(arr+shift)
    return np.abs(np.exp(log_a.mean(axis=0)) - shift)

def get_all_inst(prefix='../'):
    fpaths = glob.glob('GCB/*.col')
    names = [fpath.split('/')[-1][:-4] for fpath in fpaths]
    return names

def get_test_inst_names(method, seed, data_type, prefix='./'):
    fpaths = glob.glob(f'{prefix}/results_{data_type}/{method}/seed_{seed}/*_solving_stats.csv')
    names = [fpath.split('/')[-1][:-18] for fpath in fpaths]
    if data_type=='small':
        names = [name for name in names if name not in ["wap01a", "wap02a","wap03a", "wap04a", "C4000.5","4-FullIns_5","ash958GPIA", "C2000.5"]]
    return names

def get_test_inst_names_lp_sorted(seeds, methods, method_rets, data_type, prefix='./'):
    inst_names = get_test_inst_names(methods[0], seeds[0], data_type, prefix=prefix)
    t_sum = np.zeros(len(inst_names), dtype=float)
    for inst_id, inst_name in enumerate(inst_names):
        for seed in seeds:
            for method_ret in method_rets:
                cg_fpath = join(method_ret, f'seed_{seed}' ,f'{inst_name}_cg_stats.csv')
                with open(cg_fpath) as f:
                    t = float(f.readlines()[1].strip().split(',')[1]) 
                t_sum[inst_id] += t
    args = (-t_sum).argsort()
    inst_names = np.array(inst_names)[args].tolist()
    t_mean = t_sum / len(seeds) / len(method_rets)
    t_mean = t_mean[args]
    return inst_names, t_mean

# restricted master problem at the begining yields optimal lp bound 
# of the original master problem
def find_trivial_instances(seed, ret_path):
    fpaths = glob.glob(f'{ret_path}/*_solving_stats.csv')
    trivial_inst_names = []
    for fpath in fpaths:
        with open(fpath, 'r') as f:
            print(fpath)
            tokens = f.readlines()[1].strip().split(',')
            nb_cg_iter = int(tokens[-4])
            success = int(tokens[0] == "1")
        if success and nb_cg_iter == 1:
            trivial_inst_names.append(fpath.split('/')[-1][:-18])
    return trivial_inst_names


# no method can solve to optimality given 1000 seconds cutoff time
# def find_hard_instances():
#     inst_names = get_all_inst()
#     isHards = [False]*len(inst_names)

#     subdir_paths = [join('results', subdir_path) for subdir_path in listdir('results')]
#     for idx, inst_name in enumerate(inst_names):
        
#         optimal_count = 0
#         for subdir_path in subdir_paths:
#             solving_stats_path = join(subdir_path, f'{inst_name}_solving_stats.csv')
#             if (exists(solving_stats_path)): 
#                 with open(solving_stats_path, 'r') as f:   
#                     optimal_count += int(f.readlines()[1].strip().split(',')[0])
#             else:
#                 # print(f'solving stats file not exist: {solving_stats_path}')
#                 pass
#             isHards[idx] = optimal_count == 0

#     hard_inst_names = []
#     for idx, name in enumerate(inst_names):
#         if isHards[idx]:
#             hard_inst_names.append(name)

#     return hard_inst_names

# def trivial_insts():
#     inst_names = get_all_inst()
#     trivial_inst_names = find_trivial_instances()
#     non_trivial_inst_names = list(set(inst_names) - set(trivial_inst_names))
#     print('# inst: ', len(inst_names))
#     print('# trivial: ', len(trivial_inst_names))
#     print('# remaining: ', len(non_trivial_inst_names))
#     print(sorted(non_trivial_inst_names))