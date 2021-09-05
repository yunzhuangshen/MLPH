from utils import *
from scipy.stats import pearsonr
from scipy.stats import spearmanr

def analyse(seeds, methods, method_rets, inst_names):

    ret = [[ np.zeros(3) for j in range(len(inst_names))] for i in range(len(methods))]
    for inst_id, inst_name in enumerate(inst_names): 
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

                # log statistics
                opt = float(tot_time) < 8100 and int(optimality) == 1
                ret[method_id][inst_id][0]+=int(opt)
                ret[method_id][inst_id][1]+=float(lp_obj)
                ret[method_id][inst_id][2]+=float(tot_time)
            ret[method_id][inst_id][1:] /= len(seeds)

    diff_nb_solved_inst = np.array(ret[0])[:,0] - np.array(ret[1])[:,0]
    diff_lp = np.array(ret[0])[:,1] - np.array(ret[1])[:,1]
    tmp1 = np.array(ret[0])[:,2]
    tmp2 = np.array(ret[1])[:,2]
    print(tmp1[:10])
    tmp1[tmp1 > 1800] = 1800
    print(tmp1[:10])

    tmp2[tmp2 > 1800] = 1800
    diff_solving_time = tmp1 - tmp2

    graph_stats = []
    g_paths = [f'GCB/{g_name}.col' for g_name in inst_names]
    for g_name, g_path in zip(inst_names, g_paths):
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
        graph_stats.append([nb_node, nb_edge, density])

    graph_stats = np.array(graph_stats)
    nb_nodes, nb_edges, densities = graph_stats[:,0], graph_stats[:,1], graph_stats[:,2]

    val1, _ = pearsonr(diff_solving_time, nb_nodes)
    val2, _ = spearmanr(diff_solving_time, nb_nodes)
    print(f"correlation between solving time and nb_nodes: {val1} {val2}")
    val1, _ = pearsonr(diff_solving_time, nb_edges)
    val2, _ = spearmanr(diff_solving_time, nb_edges)
    print(f"correlation between solving time and nb_edges: {val1} {val2}")
    val1, _ = pearsonr(diff_solving_time, densities)
    val2, _ = spearmanr(diff_solving_time, densities)
    print(f"correlation between solving time and densities: {val1} {val2}")

    val1, _ = pearsonr(diff_lp, nb_nodes)
    val2, _ = spearmanr(diff_lp, nb_nodes)
    print(f"correlation between diff_lp and nb_nodes: {val1} {val2}")
    val1, _ = pearsonr(diff_lp, nb_edges)
    val2, _ = spearmanr(diff_lp, nb_edges)
    print(f"correlation between diff_lp and nb_edges: {val1} {val2}")
    val1, _ = pearsonr(diff_lp, densities)
    val2, _ = spearmanr(diff_lp, densities)
    print(f"correlation between diff_lp and densities: {val1} {val2}")


    val1, _ = pearsonr(diff_nb_solved_inst, nb_nodes)
    val2, _ = spearmanr(diff_nb_solved_inst, nb_nodes)
    print(f"correlation between diff_nb_solved_inst and nb_nodes: {val1} {val2}")
    val1, _ = pearsonr(diff_nb_solved_inst, nb_edges)
    val2, _ = spearmanr(diff_nb_solved_inst, nb_edges)
    print(f"correlation between diff_nb_solved_inst and nb_edges: {val1} {val2}")
    val1, _ = pearsonr(diff_nb_solved_inst, densities)
    val2, _ = spearmanr(diff_nb_solved_inst, densities)
    print(f"correlation between diff_nb_solved_inst and densities: {val1} {val2}")

    val2, sig = spearmanr(diff_lp, np.array(ret[0])[:,1])
    print(f"\ncorrelation between diff_lp and lp: {val2} {sig}")
if __name__ == '__main__':

    seeds = [i for i in range(1, 25)]
    methods = ['svm_cs_0', 'svm_cs_0_n']
    method_rets = [f'results_small/results_{m}' for m in methods]
    inst_names, _ = get_test_inst_names_lp_sorted(seeds, methods, method_rets, 'small')
    inst_names = inst_names[5:]

    analyse(seeds, methods, method_rets, inst_names)