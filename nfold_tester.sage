#!/usr/bin/python2
# Run as /path/to/sage nfold_tester.sage <name of pickle file> (without .pickle)
# for example /path/to/sage nfold_tester.sage QCmax_m_20_lengths_1_4_7_weigths_3_2_1_smallest_20_largest_50_slack_r_0.95
# the pickle file should contain a tuple (A, D, n, b, l, u, w, x) = inst, which is the NFoldIP instance
# Default solver is Gurobi; to use it, Sage has to be configured (https://doc.sagemath.org/html/en/thematic_tutorials/linear_programming.html#using-cplex-or-gurobi-through-sage)
# Other tested options are Coin and GLPK (GLPK is always contained in Sage, Coin has to be installed by "sage -i cbc sagelib")

# Now we set experiment parameters

# Which values of g_1 do we wish to compute for:
gc_values = [2,3,4,5,6,7,8,9,10,12,14,16,18,20,23,26,29,33,40,45,50]

# We can also compute steps of bounded \infty-norm; we call this "nginfty" method, and computing steps with bounded 1-norm we call "ng1"
methods = ["ng1"]

# Which augmentation strategies to compute for ("best" = \Gamma_best, "logarithmic" = \Gamma_2-apx, "unit" = \Gamma_unit)
gammas = ["logarithmic", "unit", "best"]

#####

import sys
import os
import subprocess
import pickle
# chdir to directory with nfoldip.pyx file
os.chdir("./")
inst_name = sys.argv[1]

#

sage.repl.load.load("nfoldip.pyx", globals())

#make new instance
pkl_file = open(inst_name+'.pickle', 'rb')
inst = pickle.load(pkl_file)
pkl_file.close()
A, D, n, b, l, u, w, x = inst


for method in methods:
    for gamma in gammas:
        for gc_val in gc_values:
            print "############################ NEW SOLVE ########################"
            name = inst_name + "_" + str(method) + "_" + gamma + "_" + str(gc_val)
            print "name:", name
            ip = NFoldIP(A, D, n, b, l, u, w, verbose = logging.INFO, graver_complexity=int(gc_val), current_solution = x,experimental=method, instancename=name, gamma=gamma, solver="Gurobi")
            ip.native_solve()

print "############################ TRUE OPTIMUM ########################"
ip.glp_solve()
