#include "training.h"
#include "MLPH.h"

namespace GCP {
    Training::Training(std::vector<std::string> training_files, std::string input_dir, double alpha, int kernel_type) :
        training_files{training_files}, input_dir{input_dir}, alpha{alpha}, kernel_type{kernel_type} {
        std::cout << "number of training graphs is: " << training_files.size() << "\n";
        construct_training_set();
    }

    void Training::construct_training_set(){
        std::string train_s = train_data_dir  + train_file_name;
        char train_data[train_s.size()+1];
        strcpy(train_data, train_s.c_str());
        
        std::ofstream train_file(train_data, std::ios::trunc);
        if (! train_file.is_open()){
            std::cout << "Cannot open the output file " <<  train_data << "\n";
            return;
        }
        train_file.close();

        std::uint64_t num0 = 0;
        std::uint64_t num1 = 0;
        for (auto d = 0u; d < training_files.size(); ++d){
            
            auto g = Instance(training_files[d], input_dir, true);
            auto n = g.size();
            
            const vector<vector<int>>& adj_list = g.adj_list;
            vector<vector<bool>> adj_matrix(n, vector<bool>(n, 0));
            for (int i = 0; i < n; ++i){
                for (int j = 0; j < adj_list[i].size(); ++j){
                    adj_matrix[i][adj_list[i][j]] = 1;
                    adj_matrix[adj_list[i][j]][i] = 1;
                }
            }

            std::cout << "number of pricing problems: " << g.nb_pp << "\n";
            for (auto inst_idx = 0; inst_idx < g.nb_pp; inst_idx++){
                std::vector<double>& cur_mis_obj_coef = g.mis_obj_coefs[inst_idx];
                std::vector<bool>& cur_mis_opt_sol = g.mis_sols[inst_idx];

                MLPH mlph(6, 0., n, n, g.degree_norm, adj_matrix, cur_mis_obj_coef, 1e8);
                mlph.random_sampling();
                mlph.compute_correlation_based_measure();
                mlph.compute_ranking_based_measure();
                train_file.open(train_data, std::ios::app);
                for (auto i = 0; i < n; ++i){
                    if (cur_mis_opt_sol[i])
                        num1++;
                    else
                        num0++;
                    float val = cur_mis_opt_sol[i];
                    float corr_norm = (mlph.max_cbm == 0) ? mlph.corr_xy[i] : mlph.corr_xy[i] / mlph.max_cbm;
                    train_file << val << " ";
                    train_file << "1:" << std::fixed << std::setprecision(6) << mlph.ranking_scores[i] / mlph.max_rbm <<  " ";
                    train_file << "2:" << std::fixed << std::setprecision(6) << corr_norm <<  " ";
                    train_file << "3:" << std::fixed << std::setprecision(6) <<  mlph.dual_values[i] / mlph.max_dual <<  " ";
                    train_file << "4:" << std::fixed << std::setprecision(6) << mlph.bound_norm[i] << " ";
                    train_file << "5:" << std::fixed << std::setprecision(6) << g.degree_norm[i] << "\n";
                }
                train_file.close();
            }
        }
        std::cout << "num0 is " << num0 << "; " << "num1 is " << num1 <<  std::endl;
        weight = alpha * num0/num1;
    }

    void Training::generate_training_model_svm(){
        std::string train_s = train_data_dir  + train_file_name;
        std::string model_s = train_data_dir + svm_train_model_name;
        char train_data[train_s.size()+1];
        char model_file[model_s.size()+1];
        strcpy(train_data, train_s.c_str());
        strcpy(model_file, model_s.c_str());
        std::cout << "weight of class 1 is " << weight << std::endl;
        std::cout << "kernel type is " << kernel_type << std::endl;
        std::cout << "output probability is " << prob << std::endl;
        svm_train_model(train_data, model_file, weight, kernel_type, prob);

        const int rem_result = remove(train_data);
        // if(rem_result == 0){
        //     std::cout << "Successfully remove training data file" << std::endl;
        // } else {
        //     std::cout << "No such training data file " << std::endl;
        // }
    }
}
