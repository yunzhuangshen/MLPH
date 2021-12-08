#include "instance.h"
#include "CG.h"
#include "training.h"
#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <random>
#include <boost/filesystem.hpp>
#include <omp.h>

using namespace GCP;
using namespace std;

static double cutoff=1800;
static double cutoff_pricer=30;
static int thread_limit=1;
string ToString(int value, int digitsCount){
    ostringstream os;
    os << setfill('0') << setw(digitsCount) << value;
    return os.str();
}

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

static double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

//extremely large instances not included: "wap03a", "wap04a", "C4000.5"
vector<string> file_name{
"1-FullIns_4", "1-FullIns_5", "2-FullIns_4", "2-FullIns_5", "3-FullIns_4", "3-FullIns_5", "4-FullIns_4", "5-FullIns_4", \
"1-Insertions_4", "1-Insertions_5", "1-Insertions_6", "2-Insertions_3", "2-Insertions_4", "2-Insertions_5", "3-Insertions_3", "3-Insertions_4", "3-Insertions_5", "4-Insertions_3", "4-Insertions_4", \
"DSJC1000.1", "DSJC1000.5", "DSJC1000.9", "DSJC125.1", "DSJC125.5", "DSJC125.9", "DSJC250.1", "DSJC250.5", "DSJC250.9", "DSJC500.1", "DSJC500.5", "DSJC500.9", \
"flat1000_50_0", "flat1000_60_0", "flat1000_76_0", "flat300_20_0", "flat300_26_0", "flat300_28_0", \
"le450_15a", "le450_15b", "le450_15c", "le450_15d", "le450_25a", "le450_25c", "le450_25d", "le450_5a", "le450_5b", "le450_5c", "le450_5d", \
"queen10_10", "queen11_11", "queen12_12", "queen13_13", "queen14_14", "queen15_15", "queen16_16", "queen8_8", "queen9_9", \
"wap01a", "wap02a", "wap05a", "wap06a", "wap07a", "wap08a", \
"mug100_1", "mug100_25", "mug88_1", "mug88_25", "myciel4", "myciel5", "myciel6", "myciel7", \
"r1000.1", "r1000.1c", "r1000.5", "r125.5", "r250.5", "will199GPIA", \
"DSJR500.1c", "DSJR500.5", "abb313GPIA", "ash608GPIA", "school1", "school1_nsh", 
};

vector<string> train_file = {"3-FullIns_4", "queen12_12", "1-Insertions_6", "mug88_25", "DSJC125.5", "flat300_20_0", "flat300_26_0", "DSJC1000.9", "DSJC250.1", "queen11_11"};

void random_ten_train_inst(){
    srand(time(NULL));
    random_shuffle(train_file.begin(), train_file.end());
    train_file.resize(10);
    cout << "training files: \n";
    for (int i = 0; i < train_file.size(); i++){
        cout << train_file[i] << " ";
    }
    cout << "\n";
}

void train_svm(){
    using namespace GCP;
    using namespace std;
    const string input_dir = "../../GCB/";
    auto training = Training(train_file, input_dir, 1, 0);
    training.generate_training_model_svm();
}

void train_logistic_model(double b0, double b1){

    using namespace GCP;
    using namespace std;
    const string input_dir = "../../GCB/";
    const int runs = 1;
    double target_sum = 0.;
    for (auto input_file_name: train_file){
        for (int i = 0; i < runs; ++i){
            const auto instance = Instance(input_file_name, input_dir, false);
            auto cg = CG(instance, 1000, 1000, 4, 1314);
            auto target = cg.optimize_LM(6, b0, b1);
            target_sum += target;
        }
    }

    string name = "../../lp_obj.txt";
    ofstream ret_file(name);
    ret_file << fixed << setprecision(5) << target_sum << "\n";
    ret_file.close();
}

int test(int method, int benchmark, int d, int column_selection, int seed, string output_dir){

    if (benchmark==0){
            cutoff=1800;
            cutoff_pricer=30;
            thread_limit=1;
    }else if(benchmark==1){
            file_name = {"wap03a", "wap04a", "C4000.5","4-FullIns_5","ash958GPIA", "C2000.5"};
            cutoff=8000;
            cutoff_pricer=150;
            thread_limit=4;
    }else{
        cout << "ERROR: unknown benchmark number\n";
        exit(-1);
    }

    for (int i = 0; i < file_name.size(); i++){
        if (d != -1 && i != d) continue;

        const string input_dir = "../../GCB/";
        string input_file_name = file_name[i];
        cout << input_file_name << endl;
        const auto instance = Instance(input_file_name, input_dir, false);
        string output_cg_filename, output_solving_filename;
        output_cg_filename = output_dir + input_file_name + "_cg_stats.csv";
        output_solving_filename = output_dir + input_file_name + "_solving_stats.csv";
        ofstream output_file_cg_stats (output_cg_filename);
        ofstream output_file_solving_stats (output_solving_filename);
        if (output_file_cg_stats.is_open()){
            output_file_cg_stats <<"ith_CG_iter,current_time,lp_obj,nnrc_cols,min_rc,mean_rc,median_rc,stdev_rc,column_selection_time,dominance_selection_time\n";
        } else{
            cout << "Cannot open the output file " + output_cg_filename << endl;
            // return 0;
        }
        if (output_file_solving_stats.is_open()){
            output_file_solving_stats << "optimality,lp_obj,tot_time,tot_cpu_time,master_duration,heur_pricing_duration,exact_pricing_duration,#CG_iter,#added_columns,#heur_success,lagrange_lower_bound" << endl;
        } else{
            cout << "Cannot open the output file " + output_solving_filename << endl;
            // return 0;
        }


        auto cg = CG(instance, cutoff, cutoff_pricer, thread_limit, seed);
        cout << "\n\n\nPROBLEM INSTANCE:" << input_file_name <<"\n";
        cout << "SOLVING ROOT LP BY CG\n"; 
        auto w0 = get_wall_time(); auto c0 = get_cpu_time();

        cg.test(method, column_selection, &output_file_cg_stats);
        
        cout << "ROOT LP SOLVED - STATS:\n";
        auto wall_time_CG = get_wall_time() - w0 - cg.time_computing_lagrangian_bound;
        auto cpu_time_CG = get_cpu_time() - c0 - cg.time_computing_lagrangian_bound;
        cout << "WALL/CPU TOTAL TIME: " << wall_time_CG << ", " << cpu_time_CG << "\n";    

        // optimality,lp_obj,
        // tot_time,tot_cpu_time, master_duration,
        // heur_pricing_duration,exact_pricing_duration,
        // #CG_iter,#added_columns,
        // #heur_success
        output_file_solving_stats << cg.lp_optimal << ","  << cg.lp_bound << ","
                << wall_time_CG << "," << cpu_time_CG << "," << cg.time_duration_master << ","
                << cg.time_duration_pricing_heur <<"," << cg.time_duration_pricing_exact << ","
                << cg.cg_iters << "," << cg.num_mis - instance.size()*10 << "," 
                << cg.num_heur_runs_success << "\n";
        
        output_file_cg_stats.close();
        output_file_solving_stats.close(); 
    }
    return 0;
}


int main(int argc, char* argv[]) {

    int mode = stoi(argv[1]);

    if (mode == 0)
        train_svm();
    else if (mode==1)
        train_logistic_model(stod(argv[2]), stod(argv[3]));
    else{
        int method = mode;
        string output_dir, benchmark, seed;
        benchmark = stoi(argv[2]) == 0 ? "small":"large";
        seed = argv[5];

        output_dir = "../results_" + benchmark + "/";
        // heuristic pricing
        int cs = 0;
        switch (method)
        {
        case 6: 
            cs = stoi(argv[4]);
            output_dir += "mlph_cs" + to_string(cs) + "/seed_" + seed + "/"; break;
        case 7: output_dir += "aco/seed_" + seed + "/"; break;
        case 8: output_dir += "tsm/seed_" + seed + "/"; break;
        case 9: output_dir += "lscc/seed_" + seed + "/"; break;
        case 10: output_dir += "fastwclq/seed_" + seed + "/"; break;
        case 11: output_dir += "gurobi/seed_" + seed + "/"; break;
        case 12: output_dir += "gurobi_heur/seed_" + seed + "/"; break;
        case 13: output_dir += "greedy/seed_" + seed + "/"; break;
        default: cout << "method is not recognized\n"; assert(false);break;
        }

        boost::filesystem::create_directories(output_dir);
        test(method, stoi(argv[2]), stoi(argv[3]), cs, stoi(seed), output_dir);
    }
    return 0;
}
