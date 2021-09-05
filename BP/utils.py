import sys, os, glob
import numpy as np
from numpy.core.fromnumeric import argsort, sort
from numpy.lib.function_base import diff
from os.path import join, isdir, isfile, exists
from os import listdir
from math import floor, ceil, isnan
from numpy import genfromtxt
from scipy.stats import ttest_ind

def get_gstats(prefix_path):
    with open(f'{prefix_path}/gprocessed.txt', 'r') as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
    
    dic = {}
    for line in lines:
        tokens = line.split(',')
        dic[tokens[0]] = (int(tokens[1]), round(float(tokens[2]), 3))
    return dic

def get_test_inst_names(method, seed, prefix_path):
    fpaths = glob.glob(f'{prefix_path}/results/{method}/{seed}/*.txt')
    names = [fpath.split('/')[-1][:-4] for fpath in fpaths]
    return names


def find_best(inst_name, seeds, methods, method_rets):
    ret = [ [0,0,[],[]] for j in range(len(methods))]
    
    for method_id, (method, method_ret) in enumerate(zip(methods, method_rets)):
        for seed in seeds:
            solving_stats_path = join(method_ret, f'{seed}', f'{inst_name}.txt')
            if not exists(solving_stats_path):
                print(f'solving stats file not exist: {solving_stats_path}')
                assert(False)
            with open(solving_stats_path, 'r') as f:  
                lines =  f.readlines()
            if len(lines) < 1:
                print(f'INCOMPLETE: {solving_stats_path}')
                assert(False)
            tokens = lines[0].strip().split(',')
            optimality, primal_bound, dual_bound, gap, nnodes, tot_time, root_time, exact_time, heur_time, fixing_col_time, primal_heur_time= tokens   
            dual_bound = ceil(float(dual_bound))
            primal_bound = float(primal_bound)
            tot_time = float(tot_time)
            optimality = int(optimality)
            hasgap=int(dual_bound != -1e20)

            ret[method_id][0]+=optimality
            ret[method_id][1]+=hasgap

            if hasgap:
                gap = (primal_bound - dual_bound)/primal_bound * 100
                ret[method_id][2].append(gap)
            if optimality:
                ret[method_id][3].append(tot_time)

    optimality_number_m1, optimality_number_m2 = ret[0][0], ret[1][0]
    hasgap_number_m1, hasgap_number_m2 = ret[0][1], ret[1][1]

    # compare optimality solved number
    if optimality_number_m1 > optimality_number_m2:
        return 0, methods[0], 'nb_solve'
    elif optimality_number_m1 < optimality_number_m2:
        return 1, methods[1], 'nb_solve'

    # compare solving time
    if optimality_number_m1 == 24 and optimality_number_m2 == 24:
        t, p = ttest_ind(ret[0][3], ret[1][3])
        if p <= 0.05:
            if t < 0:
                return 0, methods[0], 'time'
            elif t > 0:
                return 1, methods[1], 'time'

    # compare rootsolve and gap
    if hasgap_number_m1 > hasgap_number_m2:
        return 0, methods[0], 'gap'
    elif hasgap_number_m1 < hasgap_number_m2:
        return 1, methods[1], 'gap'

    t, p = ttest_ind(ret[0][2], ret[1][2])
    if p <= 0.05:
        if t < 0:
            return 0, methods[0], 'gap'
        elif t > 0:
            return 1, methods[1], 'gap'

    return -1,'None', 'None'

          
def compare(methods, method_rets, g_names, seeds, prefix_path):
    dic = get_gstats(prefix_path)

    ret=[]
    win_dic = {method:[] for method in methods}
    win_dic['none'] = []
    for g_name in g_names:

        nb_node, density = dic[g_name]
        nb_node, density = int(nb_node), float(density)
        wining_method_id, wining_method, wining_reason = find_best(g_name, seeds, methods, method_rets)
        if wining_method_id != -1:
            win_dic[wining_method].append((g_name, wining_reason))
        else:
            win_dic['none'].append(g_name)
        ret.append(f'{g_name},{nb_node},{round(density, 3)},{wining_method}')

    # with open('./graph_compare.txt', 'w') as f:
    #     f.write('x,y,label\n')
    #     for g_name, line in zip(g_names, ret):
    #         _, nb_node, density, wining_method = line.split(',')
    #         nb_node = float(nb_node)
    #         density = float(density)
    #         if wining_method != 'None':
    #             f.write(f'{nb_node},{density},{wining_method}\n')
    
    return win_dic


def find_test_insts_from_all_benchmarks(method1, method2, prefix_path):

    methods = [method1, method2]
    for method in methods:
        fpaths = glob.glob(f'{prefix_path}/results_all/{method}/1/*.txt')
        fnames = [fpath.split('/')[-1][:-4] for fpath in fpaths]

    inst_names = sorted(list(fnames))
    seeds = [1,2,3,4]
    insts_root_solved_once = []
    insts_root_not_solved = []
    insts_too_easy = []
    for inst_id, inst_name in enumerate(inst_names): 
        num_root_solved = 0; tot_optimal = 0;mean_m1=0;mean_m2=0
        for method_id, method in enumerate(methods):
            for seed in seeds:
                ret_path = f'{prefix_path}/results_all/{method}/{seed}/{inst_name}.txt'
                with open(ret_path, 'r') as f:
                    optimality, primal_bound, dual_bound, gap, nnodes, tot_time, root_time, exact_time, heur_time, fixing_col_time, primal_heur_time = f.readlines()[0].strip().split(',')
                    dual_bound = float(dual_bound)
                if dual_bound != -1e20:
                    num_root_solved+=1
                tot_optimal += int(optimality)
                if method_id==1:
                    mean_m1 += float(tot_time)
                if method_id==2:
                    mean_m2 += float(tot_time)
        
        mean_m1/=len(seeds); mean_m2/=len(seeds)
        if tot_optimal==len(seeds)*len(methods) and \
            mean_m1<10 and mean_m2 < 10:
                insts_too_easy.append(inst_name)
        else:
            if num_root_solved > 0:
                insts_root_solved_once.append(inst_name)
            else:
                insts_root_not_solved.append(inst_name)


    print("total: ", len(inst_names))
    print('too_easy: ', len(insts_too_easy))
    print('too_hard: ', len(insts_root_not_solved))
    print(insts_root_not_solved)
    print("remaining: ", len(insts_root_solved_once))
    # print(insts_root_solved_once)