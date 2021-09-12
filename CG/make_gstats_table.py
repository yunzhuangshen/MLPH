from utils import *

def make_gstats_table_all(prefix='./'):
    seeds = [i for i in range(1, 25)]
    methods = ['mlph_cs0','aco','gurobi','tsm','fastwclq','lscc']
    method_rets = [f'{prefix}/results_small/{m}' for m in methods]
    g_names, t_means = get_test_inst_names_lp_sorted(seeds, methods, method_rets, 'small', prefix=prefix)
    g_paths = [f'{prefix}/../GCB/{g_name}.col' for g_name in g_names]
    
    ret, lp_times = [], []
    for g_name, g_path, t_mean in zip(g_names, g_paths, t_means):
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

        ret.append(f'{g_name},{nb_node},{round(density, 3)},{round(t_mean,1)}')

    dest_dir = f'{prefix}/../results_cg/'
    writeto=f'{dest_dir}/gstats.txt'
    with open(writeto, 'w') as f:
        for line in ret:
            f.write(f'{line}\n')
    os.system(f'tably -n {writeto} > {dest_dir}/graph_stat_table.tex')
    os.remove(writeto)

    with open(f'{dest_dir}/graph_stat_figure_all.txt', 'w') as f:
        f.write('x,y,label\n')
        f.write('5231,0.022,large\n')
        f.write('4730,0.026,large\n')
        f.write('4146,0.009,large\n')
        f.write('4000,0.500,large\n')
        f.write('2464,0.037,large\n')
        f.write('2368,0.04,large\n')
        f.write('2000,0.500,large\n')
        f.write('1916,0.007,large\n')

        for line in ret:
            _, nb_node, density, lp_time = line.split(',')
            nb_node = float(nb_node)
            density = float(density)
            lp_time = float(lp_time)
            label='small'
            f.write(f'{nb_node},{density},{label}\n')
            


def make_gstats_table_small(prefix='./'):
    seeds = [i for i in range(1, 25)]
    methods = ['mlph_cs0','aco','gurobi','tsm','fastwclq','lscc']
    method_rets = [f'{prefix}/results_small/{m}' for m in methods]
    g_names, t_means = get_test_inst_names_lp_sorted(seeds, methods, method_rets, 'small', prefix=prefix)
    g_paths = [f'{prefix}/../GCB/{g_name}.col' for g_name in g_names]
    
    ret, lp_times = [], []
    for g_name, g_path, t_mean in zip(g_names, g_paths, t_means):
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

        ret.append(f'{g_name},{nb_node},{round(density, 3)},{round(t_mean,1)}')
    
    dest_dir = f'{prefix}/../results_cg/'
    with open(f'{dest_dir}/graph_stat_figure_small.txt', 'w') as f:
        f.write('x,y,label\n')
        for line in ret:
            _, nb_node, density, lp_time = line.split(',')
            nb_node = float(nb_node)
            density = float(density)
            f.write(f'{nb_node},{density}\n')

if __name__ == '__main__':
    make_gstats_table_small()