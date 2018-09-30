import os
import pickle
from collections import defaultdict
from itertools import product
import pandas
import seaborn as sns
import numpy as np


def load_data(LOGDIR, stripped=True):
    d = {}
    if stripped:
        condition = ".stripped.pickle"
    else:
        condition = ".pickle"
    for root, dirs, files in os.walk(LOGDIR):
        for name in files:
            if condition in name:
                p=open(os.path.join(root, name), "rb")
                try:
                    log = pickle.load(p)
                except:
                    print("failed:",name)
                p.close()
                log["directory"] = root
                d[log["instancename"]] = log
    return d

def master_df(d):
    # First precompute for each directory the minimum number
    # of iterations needed to reach optimum
    # (called it_min in the docs)
    dirs = defaultdict(lambda: float("Inf"))
    for inst in d:
        if d[inst]["native_solution"] == d[inst]["glpk_solution"]:
            iterations = sum(sum((1 for gamma in d[inst]["iterations"][it]["gammas"])) for it in d[inst]["iterations"])
            dirs[d[inst]["directory"]] = min(dirs[d[inst]["directory"]], iterations)
    
    
    breakdown = {"dimension": [], "Delta": [], "time": [], "time type": [], 
                 "gc": [], "iter": [], "obj type": [], "obj": [],
                 "Gamma": [], "gamma": [], "last_gamma": [],
                "n": [], "t": [], "r": [],
                "directory": [], "gap": [], "convergence": [], "fail": []}
    for inst in d:
        it_index = 0
        
        total = d[inst]["end_time"] - d[inst]["start_time"]
        augip_init = d[inst]["augip_init_time"]
        gurobi_solve = d[inst]["glpk_solve_time"]*10**2
        gurobi_construct_solve = d[inst]["glpk_construct_and_solve"]
        augip_solve = 0.0
        augip_solves = 0
        Delta, Gamma = d[inst]["Delta"], d[inst]["gamma"]
        n, r, t, directory = d[inst]["inst"]["n"], d[inst]["inst"]["r"], d[inst]["inst"]["t"], d[inst]["directory"]
        gc = d[inst]["graver_complexity"]
        aug_opt = d[inst]["native_solution"]
        gurobi_opt = d[inst]["glpk_solution"]
        gurobi_fail = (gurobi_opt is False)
        
        if aug_opt > gurobi_opt:
            it_g = float("Inf")
        else:
            it_g = sum(sum((1 for gamma in d[inst]["iterations"][it]["gammas"])) for it in d[inst]["iterations"])
        it_min = dirs[d[inst]["directory"]]
        convergence = it_min / float(it_g)
        
        for iteration in d[inst]["iterations"]:
            it = d[inst]["iterations"][iteration]
            iter_opt = min(d[inst]["iterations"][iteration]["gammas"][ga]["obj"] for ga in d[inst]["iterations"][iteration]["gammas"])
            for gamma in it["gammas"]:
                it_index += 1
                augip_solve += it["gammas"][gamma]["solve_time"]
                augip_solves += 1
                augip_solve_average = (augip_solve/augip_solves)*10**2
                val = it["gammas"][gamma]["obj"]
                last_gamma = (gamma == max(it["gammas"].keys()))
                aug_fail = (it["gammas"][gamma]["solve_success"] is False)
                for (value, typ) in [(val, "actual"), (iter_opt, "min")]:
                    breakdown["obj"].append(value)
                    breakdown["obj type"].append(typ)
                    breakdown["dimension"].append(d[inst]["dimension"])
                    breakdown["Delta"].append(d[inst]["Delta"])
                    breakdown["time"].append(it["gammas"][gamma]["solve_time"]*10**2)
                    breakdown["time type"].append("aug scaled 10^2")
                    breakdown["Gamma"].append(Gamma)
                    breakdown["gamma"].append(gamma)
                    breakdown["gc"].append(gc)
                    breakdown["iter"].append(it_index)
                    breakdown["last_gamma"].append(last_gamma)
                    breakdown["n"].append(n)
                    breakdown["r"].append(r)
                    breakdown["t"].append(t)
                    breakdown["directory"].append(directory)
                    breakdown["gap"].append(value - gurobi_opt)
                    breakdown["convergence"].append(convergence)
                    breakdown["fail"].append(aug_fail)

        
        if augip_solves == 0:
            continue
        augip_solve_average = (augip_solve/augip_solves)*10**2

        for (cat,time_val, obj_val, obj_typ, fail) in [("total", total, aug_opt, "aug", False),
                                                       ("aug total", augip_solve, aug_opt, "aug", False),
                        ("aug init", augip_init, aug_opt, "aug", False), ("aug average scaled 10^2", augip_solve_average, aug_opt, "aug", False),
                         ("gurobi solve scaled 10^2", gurobi_solve, gurobi_opt, "gurobi", gurobi_fail),
                          ("gurobi construct & solve", gurobi_construct_solve, gurobi_opt, "gurobi", gurobi_fail)]:
            breakdown["dimension"].append(d[inst]["dimension"])
            breakdown["Delta"].append(Delta)
            breakdown["time"].append(time_val)
            breakdown["time type"].append(cat)
            breakdown["Gamma"].append(Gamma)
            breakdown["gamma"].append(-1)
            breakdown["gc"].append(gc)
            breakdown["obj"].append(obj_val)
            breakdown["obj type"].append(obj_typ)
            breakdown["iter"].append(-1)
            breakdown["last_gamma"].append(False)
            breakdown["n"].append(n)
            breakdown["r"].append(r)
            breakdown["t"].append(t)
            breakdown["directory"].append(directory)
            breakdown["gap"].append(obj_val - gurobi_opt)
            breakdown["convergence"].append(convergence)
            breakdown["fail"].append(fail)
            
    data = pandas.DataFrame(breakdown)
    return data
        
#def heatmap_data(df, row, col, val):
    ## generates data for heatmap indexed by row and col,
    ## each cell of the heatmap is mean of val
    ## over all points from that row-col pair.
    #base = defaultdict(list)
    #rows = list(set(df[row]))
    #cols = list(set(df[col]))
    #rows.sort()
    #cols.sort()
    #for (r,c) in product(rows, cols):
        #base[(r,c)] = df[df[row] == r][df[col] == c][val].agg(["mean"])[0]

    #table = []
    #for r in rows:
        ## this will be one row of the table
        #table_row = []
        #for c in cols:
            #if (r,c,) in base:
                #table_row.append(base[(r,c)])
        #table.append(table_row)

    #data = pandas.DataFrame(table, index=rows, columns=cols)
    #return data

def heatmap_data(df, row, col, val):
    return pandas.pivot_table(df, values=val, index=row, columns=col, aggfunc=np.mean)
    
        

