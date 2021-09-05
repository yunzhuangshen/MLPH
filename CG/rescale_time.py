import sys, os
from os import fpathconf, listdir
from os.path import isfile, join
import numpy as np

if __name__ == '__main__':
    # input_dirs = ['small/lp-curve', 'small/solving-curve']; divideby=100

    input_dirs = ['large/lp-curve', 'large/solving-curve',\
                    'large-cs/lp-curve', 'large-cs/solving-curve']; divideby=1000

    input_dirs = [f'data/{input_dir}' for input_dir in input_dirs]

    for input_dir in input_dirs:
        fpaths = [join(input_dir, f) for f in listdir(input_dir) if isfile(join(input_dir, f))]
        for fpath in fpaths:
            M = np.loadtxt(fpath,delimiter=",",skiprows=1,dtype=float)
            M[:,0]/=divideby
            np.savetxt(fpath, M, delimiter=",")
            with open(fpath, 'r') as f:
                lines = f.readlines()
                lines.insert(0,'x,y\n')
            with open(fpath,'w') as f:
                for line in lines:
                    f.write(line)