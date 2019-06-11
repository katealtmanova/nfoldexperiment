#!/aux/alquaknaa/sage/SageMath/sage
##!/usr/bin/python2

import argparse

#list_machines = range(10,50,10) #+[125,150,175,200]+[250,300]
#count_for_each_p = 4
#number_of_job_types = 4
#p_s = [5,6,7,8,9,10]
#slacks = [0.5, 0.6, 0.65, 0.7, 0.75, 0.8]

parser = argparse.ArgumentParser(description='n-fold IP scheduling tester')
parser.add_argument('--instance_type', nargs='?', default='sched', type=str, help='sched or cs, determines whether we generate instances of Q||C_max or Closest String.')
parser.add_argument('--logdir', nargs='?', default="logs", type=str, help="Directory into which logs will be saved. Useful for separating different projects.")
parser.add_argument('--gc_values', nargs='*', default=[4,8,12,20,30,40,50,75,100], type=int, help='Values of tuning parameter.')
parser.add_argument('--gammas', nargs="*", default=["log2"], type=str, help="Gammas to test. Options are log2, log5, log10, unit, best")

parser.add_argument('--disable_nfold', dest='nfold', action='store_const',
                    const=False, default=True,
                    help='Disable the iterative solver and only run Gurobi on the generated instances.')
parser.add_argument('--augip_timelimit', nargs="?", default=60, type=int, help="Timelimit for AugIP")
parser.add_argument('--milp_timelimit', nargs="?", default=60, type=int, help="Timelimit for Gurobi solve of the full instance")


# Specific for sched
parser.add_argument('--count_for_each_p', nargs="?", default=3, type=int, help="How many samples for each Delta (sched).")
parser.add_argument('--p_s', nargs='*', default=[5,6,7,8,9,10,11,12,13], type=int, help='Pick processing times from the set of first P primes (sched).')
parser.add_argument('--slacks', nargs='*', default=[0.6, 0.7, 0.8], type=float, help='Slack of instance (sched).')
parser.add_argument('--machines', nargs='*', default=[10,20,30,40,50,60,70,80,90,100], type=int, help='Numebr of machines (sched).')
parser.add_argument('--number_job_types', nargs="*", default=[4], type=int, help="Number of job types aka $t$ in the n-fold formulation (sched).")


# Specific for CS
parser.add_argument('--str_len', nargs="*", default=[500,1000,2000,4000,8000,16000], type=int, help="Length of generated strings (cs).")
parser.add_argument('--str_num', nargs='*', default=[3,4,5,6], type=int, help='Numebr of strings (cs).')
parser.add_argument('--ratio', nargs='*', default=[2,3,4,7,10,15], type=float, help='Ratio n/alpha, so the random instance has alpha changes compared with the generated master string (cs).')
parser.add_argument('--sigma', nargs='*', default=[2,3,4,5], type=int, help='Size of alphabet (cs).')
parser.add_argument('--distance_factor', nargs="*", default=[0.1,0.15,0.2,0.25,0.3,0.5,0.7], type=float, help="Distance factor, float between 0 and 1. If 0 the instance is almost surely infeasible, if 1 it should be always feasible.")

args = parser.parse_args()


# Now we set experiment parameters

# Which values of g_1 do we wish to compute for:
#gc_values = [2,3,4,5,6,7,8,9,10,12,14,16,18,20,23,26,29,33,40,45,50]
#gc_values = [4,8,12,20,30,40,50,75,100]

# We can also compute steps of bounded \infty-norm; we call this "nginfty" method, and computing steps with bounded 1-norm we call "ng1"
methods = ["ng1"]

# Which augmentation strategies to compute for ("best" = \Gamma_best, "logarithmic" = \Gamma_2-apx, "unit" = \Gamma_unit)
# gammas = ["logarithmic", "unit", "best"]




#####

import sys
import os
import subprocess
import pickle
import random
import copy
import hashlib
import itertools
from sage.sets.primes import Primes
# chdir to directory with nfoldip.pyx file

from sched_instance_generator import gen_sched_instance
from cs_instance_generator import gen_cs_instance




def gen_sched_instances():
    P = Primes()

    instances = []

    primes_up_to = [P.unrank(i) for i in range(max(args.p_s)+1)]
    for p in args.p_s:
        Delta = int(P.unrank(p))
        for t in args.number_job_types:
            for i in range(args.count_for_each_p):
                for sl in args.slacks:
                    inst = {}
                    job_ls = random.sample(primes_up_to[:p], t-1) + [Delta] # we always add the largest prime so that all samples in this batch have the same Delta
                    job_ls.sort()
                    inst["job_lengths"] = copy.deepcopy(job_ls)
                    job_ls.reverse()
                    inst["job_weights"] = job_ls
                    inst["smallest"] = 5*Delta
                    inst["largest"] = int(Delta**(2.5))
                    inst["obj"] = "total_jobs"
                    inst["slack_r"] = sl
                    inst["number_job_types"] = t
                    for m in args.machines:
                        inst_m = copy.deepcopy(inst)
                        inst_m["m"] = m
                        instances.append(inst_m)
    def m_plus_t(inst):
        return inst["m"] + inst["number_job_types"]
    instances.sort(key=m_plus_t)
    return instances


def gen_cs_instances():
    instances = []
    for (n,k,r,s,dist_factor) in itertools.product(args.str_len, args.str_num,
                                                   args.ratio, args.sigma, args.distance_factor):
        inst = dict(zip(("n","k","r","sigma","dist_factor"),(n,k,r,s,dist_factor)))
        instances.append(inst)
    return instances

if args.instance_type == "sched":
    sched_instances = gen_sched_instances()
elif args.instance_type == "cs":
    cs_instances = gen_cs_instances()

def sage_int_matrix_to_list_of_tuples(M):
    # Converts a sage integer matrix to a list of rows,
    # each row being a tuple of ints
    # thus matrix(rs) would return M
    rs = M.rows()
    rs= list(tuple((int(ri) for ri in r)) for r in rs)
    return rs



sage.repl.load.load("nfoldip.pyx", globals())

def run_tests(inst):
    method = "ng1"
    A, D, n, b, l, u, w, x, inst_name = inst
    print("u",u)
    print("x",x)
    inst_name += "_" + str(random.randint(10000,99999))
    Delta = int(max((A.numpy().max(), D.numpy().max())))
    dimension = n*D.ncols() # n*t
    path = os.path.join("./", args.logdir, str(dimension), str(Delta), inst_name)
    os.makedirs(path)
    os.chdir(path)
    if args.nfold:
        for gc_val in args.gc_values:
            for gamma in args.gammas:
                print "############################ NEW SOLVE ########################"
                name = inst_name + "_" + str(gc_val) + "_" + gamma
                name = str(hashlib.sha224(name).hexdigest()) # Using this because we ran into names too long
                print "name:", name
                ip = NFoldIP(A, D, n, b, l, u, w, verbose = logging.INFO, graver_complexity=int(gc_val), current_solution = x,
                             experimental=method, instancename=name, gamma=gamma, solver="Gurobi", augip_timelimit=args.augip_timelimit, milp_timelimit=args.milp_timelimit)
                ip.native_solve()
                ip.glpk_solve()
                ip.fulllog["nfold"] = args.nfold
                fulllog_f = open(name+".pickle", "w")
                pickle.dump(ip.fulllog, fulllog_f)
                fulllog_f.close()
    else:
        print "############################ NEW SOLVE ########################"
        name = inst_name + "_onlymilp"
        name = str(hashlib.sha224(name).hexdigest())
        print "name:", name
        ip = NFoldIP(A, D, n, b, l, u, w, verbose = logging.INFO, graver_complexity=int(1), current_solution = x,
                     experimental=method, instancename=name, solver="Gurobi", augip_timelimit=args.augip_timelimit, milp_timelimit=args.milp_timelimit)
        ip.glpk_solve()
        ip.fulllog["nfold"] = args.nfold
        fulllog_f = open(name+".pickle", "w")
        pickle.dump(ip.fulllog, fulllog_f)
        fulllog_f.close()


    os.chdir("../../../../")


if args.instance_type == "sched":
    for (i,inst) in enumerate(sched_instances): # -> FIXME instances when we are ready to go
        print("####### Testing instance",i,"/",len(sched_instances))
        ip_instance = gen_sched_instance(inst["job_lengths"], inst["job_weights"], inst["smallest"], inst["largest"], inst["m"], slack_r=inst["slack_r"], obj=inst["obj"])
        run_tests(ip_instance)
elif args.instance_type == "cs":
    for inst in cs_instances:
        df = "{0:.2f}".format(float(inst["dist_factor"]))
        inst_name = "_".join(("cs", str(inst["n"]), str(inst["k"]), str(inst["r"]), str(inst["sigma"]), df))
        ip_instance = gen_cs_instance(inst["n"], inst["k"], inst["r"], inst["sigma"], inst["dist_factor"]) + (inst_name,)
        
        run_tests(ip_instance)
    





