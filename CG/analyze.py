from utils import *
from make_compare_table import make_compare_table
from make_solving_table import make_solving_table
from make_rc_table import make_rc_table
from make_cg_table import make_cg_table
import pathlib

def solving_curve(seeds, methods, method_rets, inst_names, outdir_name, interval=5):

    plots = {}
    if 'small' in outdir_name:
        tot_time = 1800
    else:
        tot_time = 8000
    nslots = tot_time//interval
    for method in methods:
        plots[method] = [0 for i in range(nslots)] # time slot every 5 seconds for 1800 seconds cutoff
            
    for method, method_ret_fname in zip(methods, method_rets):
        for inst_name in inst_names:
            for seed in seeds:
                subdir_path = join(method_ret_fname, f'seed_{seed}')
                solving_stats_path = join(subdir_path, f'{inst_name}_solving_stats.csv')

                if not exists(solving_stats_path):
                    print(subdir_path)
                    print(f'solving stats file not exist: {solving_stats_path}')
                    assert(False)
                with open(solving_stats_path, 'r') as f:  
                    lines =  f.readlines()
                if len(lines) < 2:
                    print(subdir_path)
                    print(f'ERROR at solving stats file: {solving_stats_path}')
                    assert(False)
                tokens = lines[1].strip().split(',')
                optimal = tokens[0]=='1'
                if optimal:
                    tokens = lines[1].strip().split(',')
                    idx = floor(float(tokens[2]) // interval)
                    for i in range(idx, nslots):
                        plots[method][i]+=1
    
    os.makedirs(f'{outdir_name}/solving-curve',exist_ok=True)
    for method, curve in plots.items():
        xs = [interval*i for i in range(1, nslots+1)]
        ys = curve
        assert(len(xs) == len(ys))
        with open(f'{outdir_name}/solving-curve/{method.replace("_","-")}.txt', 'w') as f:
            f.write('x,y\n')
            for i in range(len(xs)):
                f.write(f'{xs[i]},{ys[i]}\n')


def solving_curve_cg(seeds, methods, method_rets, inst_names, outdir_name, ncg_iters=10):

    plots = {}
    ncg_iters = 5
    for method in methods:
        plots[method] = [0 for i in range(ncg_iters)] # time slot every 5 seconds for 1800 seconds cutoff
            
    for method, method_ret_fname in zip(methods, method_rets):
        for inst_name in inst_names:
            for seed in seeds:
                subdir_path = join(method_ret_fname, f'seed_{seed}')
                solving_stats_path = join(subdir_path, f'{inst_name}_solving_stats.csv')
                cg_stats_path = join(subdir_path, f'{inst_name}_cg_stats.csv')

                if not exists(solving_stats_path):
                    print(subdir_path)
                    print(f'solving stats file not exist: {solving_stats_path}')
                    assert(False)
                with open(solving_stats_path, 'r') as f:  
                    lines =  f.readlines()
                if len(lines) < 2:
                    print(subdir_path)
                    print(f'ERROR at solving stats file: {solving_stats_path}')
                    assert(False)
                tokens = lines[1].strip().split(',')
                optimal = tokens[0]=='1'
                if optimal:
                    with open(cg_stats_path, 'r') as f:  
                        solved_ith_cg_iteration =  len(f.readlines()[1:])-1-1

                    for i in range(solved_ith_cg_iteration, ncg_iters):
                        plots[method][i]+=1
    
    os.makedirs(f'{outdir_name}/solving-curve-cg',exist_ok=True)
    for method, curve in plots.items():
        xs = [i for i in range(1, ncg_iters+1)]
        ys = curve
        assert(len(xs) == len(ys))
        with open(f'{outdir_name}/solving-curve-cg/{method.replace("_","-")}.txt', 'w') as f:
            f.write('x,y\n')
            for i in range(len(xs)):
                f.write(f'{xs[i]},{ys[i]}\n')
                
ARITHMATIC=0;GEOMETRIC=1
def lp_curve(seeds, methods, method_rets, inst_names, mean_type, outdir_name, interval=5):

    stats_dict = {}
    if 'small' in outdir_name:
        tot_time = 1800
    else:
        tot_time = 8000
    nslots = tot_time//interval
    for method in methods:
        stats_dict[method] = [[] for i in range(nslots)] # time slot every 5 seconds for 1800 seconds cutoff

    for method, method_ret_fname in zip(methods, method_rets):
        for inst_name in inst_names:
            for seed in seeds:
                subdir_path = join(method_ret_fname, f'seed_{seed}')
                cg_stats_path = join(subdir_path, f'{inst_name}_cg_stats.csv')

                if not exists(cg_stats_path):
                    print(f'solving stats file not exist: {cg_stats_path}')
                    assert(False)
                with open(cg_stats_path, 'r') as f:  
                    lines =  f.readlines()[1:]
                
                indices = [floor(float(line.strip().split(',')[1]) //interval) for line in lines]
                lps = [float(line.strip().split(',')[2]) for line in lines]

                indices.insert(0, 0); lps.insert(0,lps[0])

                # remove data out of time limit
                tmp_indices=[]; tmp_lps=[]
                for i in range(len(indices)):
                    if indices[i] >= nslots:
                        break
                    tmp_indices.append(indices[i])
                    tmp_lps.append(lps[i])

                indices, lps = tmp_indices, tmp_lps
                assert(len(indices) == len(lps))
                # remove duplicates
                tmp_indices=[]; tmp_lps=[]
                for i in range(len(indices)-1):
                    if indices[i] < indices[i+1]:
                        tmp_indices.append(indices[i])
                        tmp_lps.append(lps[i])

                tmp_indices.append(indices[-1])
                tmp_lps.append(lps[-1])
                if tmp_indices[-1]!=nslots-1:
                    tmp_indices.append(nslots-1)
                    tmp_lps.append(lps[-1])

                indices, lps = tmp_indices, tmp_lps
                for i in range(len(indices)-1):
                    _from = indices[i]
                    _to = indices[i+1]
                    # print(_from, _to, len(stats_dict[method]))
                    for j in range(_from, _to):
                        stats_dict[method][j].append(lps[i])
                stats_dict[method][-1].append(lps[-1])


    os.makedirs(f'{outdir_name}/lp-curve',exist_ok=True)
    for method, curve in stats_dict.items():
        xs = [interval*i for i in range(1, nslots+1)]
        ys = []
        for arr in curve:
            assert(len(arr)==len(seeds)*len(inst_names))
            if mean_type == ARITHMATIC:
                ys.append(np.mean(arr))
            else:
                ys.append(gmean(arr))
        
        with open(f'{outdir_name}/lp-curve/{method.replace("_","-")}.txt', 'w') as f:
            f.write('x,y\n')
            for i in range(len(xs)):
                f.write(f'{xs[i]},{round(ys[i],5)}\n')

          
def lp_curve_cg(seeds, methods, method_rets, inst_names, mean_type, outdir_name, ncg_iters=10):

    stats_dict = {}
    for method in methods:
        stats_dict[method] = [[] for i in range(ncg_iters)] # time slot every 5 seconds for 1800 seconds cutoff

    for method, method_ret_fname in zip(methods, method_rets):
        for inst_name in inst_names:
            for seed in seeds:
                subdir_path = join(method_ret_fname, f'seed_{seed}')
                cg_stats_path = join(subdir_path, f'{inst_name}_cg_stats.csv')

                if not exists(cg_stats_path):
                    print(f'solving stats file not exist: {cg_stats_path}')
                    assert(False)
                with open(cg_stats_path, 'r') as f:  
                    lines =  f.readlines()[1:]
                
                lps = [float(line.strip().split(',')[2]) for line in lines]
                if (len(lps)==0):
                    print(cg_stats_path)
                if ncg_iters < len(lps):
                    lps=lps[:ncg_iters]
                else:
                    while ncg_iters > len(lps):
                        lps.append(lps[-1]);                        
                    assert(ncg_iters==len(lps))

                for i in range(ncg_iters):
                    stats_dict[method][i].append(lps[i])

    os.makedirs(f'{outdir_name}/lp-curve-cg',exist_ok=True)
    for method, curve in stats_dict.items():
        xs = [i for i in range(1, ncg_iters+1)]
        ys = []
        for arr in curve:
            assert(len(arr)==len(seeds)*len(inst_names))

            if mean_type == ARITHMATIC:
                ys.append(np.mean(arr))
            else:
                ys.append(gmean(arr))
        
        with open(f'{outdir_name}/lp-curve-cg/{method.replace("_","-")}.txt', 'w') as f:
            f.write('x,y\n')
            for i in range(len(xs)):
                f.write(f'{xs[i]},{round(ys[i],5)}\n')

def analyse_small(prefix):

    seeds = [i for i in range(1, 25)]
    methods = ['mlph_cs0','aco','gurobi','gurobi_heur','tsm','fastwclq','lscc']
    method_rets = [f'{prefix}/results_small/{m}' for m in methods]
    inst_names, t_means = get_test_inst_names_lp_sorted(seeds, methods, method_rets, 'small', prefix=prefix)
    dest_dir = f'{prefix}/../results_cg/small/'
    os.makedirs(dest_dir,exist_ok=True)
    make_solving_table(seeds, methods, method_rets, inst_names, dest_dir)
    make_rc_table(seeds, methods, method_rets, inst_names, dest_dir)
    make_cg_table(seeds, methods, method_rets, inst_names, dest_dir)
    solving_curve(seeds, methods, method_rets, inst_names, dest_dir)
    lp_curve(seeds, methods, method_rets, inst_names, GEOMETRIC, dest_dir)

def analyse_large(prefix):
    seeds = [i for i in range(1, 25)]
    inst_names = get_test_inst_names( 'mlph_cs0', seeds[0], 'large', prefix=prefix)
    methods = ['mlph_cs0', 'aco', 'gurobi', 'gurobi_heur', 'fastwclq']
    method_rets = [f'{prefix}/results_large/{m}' for m in methods]
    dest_dir = f'{prefix}/../results_cg/large/'
    os.makedirs(dest_dir,exist_ok=True)
    make_solving_table(seeds, methods, method_rets, inst_names, dest_dir)
    make_rc_table(seeds, methods, method_rets, inst_names, dest_dir)
    make_cg_table(seeds, methods, method_rets, inst_names, dest_dir)
    solving_curve(seeds, methods, method_rets, inst_names, dest_dir)
    lp_curve(seeds, methods, method_rets, inst_names, GEOMETRIC, dest_dir)

def analyse_cs_large(prefix):
    seeds = [i for i in range(1, 9)]
    inst_names = get_test_inst_names( 'mlph_cs0', seeds[0], 'large', prefix=prefix)
    
    # compare column selection method
    methods=[]
    for i in [0,1,5]:
        methods.append(f'mlph_cs{i}')
    method_rets = [f'{prefix}/results_large/{m}' for m in methods]
    dest_dir = f'{prefix}/../results_cg/cs-large'
    os.makedirs(dest_dir,exist_ok=True)

    # solving_curve(seeds, methods, method_rets, inst_names, dest_dir)
    # solving_curve_cg(seeds, methods, method_rets, inst_names, dest_dir)
    lp_curve(seeds, methods, method_rets, inst_names, GEOMETRIC, dest_dir)
    lp_curve_cg(seeds, methods, method_rets, inst_names, GEOMETRIC, dest_dir)

def analyze(prefix='./'):
    analyse_small(prefix)
    analyse_large(prefix)
    analyse_cs_large(prefix)
    make_compare_table(prefix=prefix)
if __name__ == '__main__':
    analyze()