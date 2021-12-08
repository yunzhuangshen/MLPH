import sys, os
import numpy as np
from bayes_opt import BayesianOptimization
import sys,os
from subprocess import check_call
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
from bayes_opt.util import load_logs
import datetime
from time import time



def train_svm_linear():
    def dual_to_primal():
        dual_path = 'CG/train_data/svm_train_model'
        primal_path = 'svm.param'

        with open(dual_path, 'r') as f:
            lines = f.readlines()
            lines = [line.strip() for line in lines if len(line)!=0]
        
        nsv = int(lines[3].split(' ')[1]) 
        b = -float(lines[4].split(' ')[1])
        svs = lines[8:]
        assert(len(svs) == nsv)
        w = np.zeros(5, dtype=float)
        for sv in svs:
            tokens = sv.split(' ')
            coef = float(tokens[0])
            for i, feature in enumerate(tokens[1:]):
                feature = float(feature.split(':')[1])
                w[i] += feature*coef
        
        with open(primal_path, 'w+') as f:
            for weight in w.tolist():
                f.write(f'{weight}\n')
            f.write(f"{b}\n")

    os.system(f'cd ./CG/build/ && ./CG 0 > tmp.txt')
    dual_to_primal()



def optimize_logistic_model():
    def black_box_function(b0, b1):
        """Function with unknown internals we wish to maximize."""
        b0 = str(b0.item())
        b1 = str(b1.item())
        t0 = time()
        os.system(f'cd ./CG/build/ && ./CG 1 {b0} {b1}')
        print(f'time used current run: {t0-time()}')
        with open(f'./lp_obj.txt', 'r') as f:
            obj = float(f.readlines()[0].strip())
        return -obj

    # Bounded region of parameter space
    pbounds = {
                'b0': (-15, 0.0), 
                'b1': (-30, 0),
                }

    optimizer = BayesianOptimization(
        f=black_box_function,
        pbounds=pbounds,
        random_state=1314,
    )

    logger = JSONLogger(path=f"BO.json")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)
    optimizer.maximize(
        init_points=150,
        n_iter=150,
    )
    if os.path.exists('./lp_obj.txt'):
        os.remove('./lp_obj.txt')

    b0, b1 = optimizer.max['params']['b0'], optimizer.max['params']['b1']
    with open("lm.param", 'w+') as f:
        f.write(f'{b0}\n{b1}\n')
    


if __name__ == '__main__':

    os.system(f'mkdir CG/build; cd CG/build && cmake ../ && make')
    train_svm_linear()
    optimize_logistic_model()