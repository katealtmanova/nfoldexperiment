from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.rings.infinity import Infinity
from sage.matrix.special import block_matrix
from sage.matrix.special import diagonal_matrix
from sage.matrix.special import zero_matrix
from sage.matrix.special import identity_matrix
from sage.functions.other import floor, ceil
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.structure.sage_object import SageObject
from sage.interfaces.four_ti_2 import four_ti_2
from sage.numerical.mip import MIPSolverException

from itertools import combinations_with_replacement
from math import log

import random, collections
def gen_sched_instance(job_lengths, job_weights, smallest, largest, m, slack_r=0.9, dist=None, obj="total_jobs"):
    # job_lengths is list of ints
    # job_weights is a list of ints determining the distribution of jobs
    # smallest, largest are ints - Cmax for smallest and longest machine
    # m is number of machines
    # slack is the total slack of the instance (total capacities of machines - total sizes of jobs);
    # slack_r is the RATIO of total capacity / total sizes
    # obj is ["total_length", "total_jobs", "order"]
    
    inst_name = "QCmax_m_" + str(m)
    inst_name += "_lengths_" + "_".join((str(j) for j in job_lengths))
    inst_name += "_weigths_" + "_".join((str(j) for j in job_weights))
    inst_name += "_smallest_" + str(smallest) + "_largest_" + str(largest)
    inst_name += "_slack_r_" + "%.2f" % slack_r
    inst_name += "_obj_" + obj
    
    print inst_name
    
    # Algorithm:
    # 1. pull m machine capacities from some distribution (between smallest, largest)
    # 2. compute total capacity
    # 3. compute total sizes
    # 4. keep randomly pulling job lengths until we reach the slack; once it would be violated, terminate.
    
    # Default distribution is uniform.
    M = [random.randrange(smallest, largest+1) for i in range(m)]
    M.sort()
    
    total_cap = sum(M)
    total_size_limit = int(total_cap * slack_r)
    
    jobs = collections.defaultdict(int)
    
    total_size = 0
    weighted_job_lengths = []
    for (i,j) in enumerate(job_lengths):
        for k in range(job_weights[i]):
            weighted_job_lengths.append(j)
    
    
    #weighted_job_lengths = [j for k in range(job_weights[i]) for (i,j) in enumerate(job_lengths)]
    
    while True:
        j = random.choice(weighted_job_lengths)
        total_size += j
        if total_size > total_size_limit:
            break
        else:
            jobs[j] += 1
    job_vector = tuple()
    for l in job_lengths:
        job_vector += (l,)
            
    ##### NOW generate n-fold instance        
    
    r = len(job_lengths)
    t = r + 1
    A = job_lengths + [1]
    D = []
    for i in range(r):
        row = [0]*i + [1] + [0]*(r-i)
        D.append(row)
        
    l = [(0,)*t]*(m+1)
    u = [job_vector + (total_size,)] + [job_vector + (M[i],) for i in range(m)]
    if obj == "total_jobs":
        w = [(1,)*r+(0,)] + [(0,)*t]*m
    b = [job_vector] + [(total_size,)] + [(M[i],) for i in range(len(M))]
    x = [job_vector + (0,)] + [(0,)*r + (M[i],) for i in range(len(M))]
    
    l = [vector(ll) for ll in l]
    u = [vector(uu) for uu in u]
    w = [vector(ww) for ww in w]
    x = [vector(xx) for xx in x]
    
    return matrix(A), matrix(D), m+1, b, l, u, w, x, inst_name
