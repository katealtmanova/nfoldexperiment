# Evaluating and Tuning n-Fold Integer Programming
Supplementary material.

## Code

- The main class is [nfoldip.pyx](./nfoldip.pyx)
- Its documentation is in [docs_NFoldIP.md](./docs_NFoldIP.md)
- A tester script, which takes in a pickle file with an n-fold IP instance and runs the solver for given values of the parameter g1, augmentation strategies etc. [nfold_tester.sage](./nfold_tester.sage)
- A script for generation of random Closest String instances [./jupyter/cs_random_instance_generator.ipynb](./jupyter/cs_random_instance_generator.ipynb) and functions for generating n-fold IP instances from the standard format [./sage/csp_instance_generation.sws](./sage/csp_instance_generation.sws)
- A script for generation of random Scheduling instances [./sage/random_scheduling_instances.sws](./sage/random_scheduling_instances.sws)
- [Outer loop plotter](./jupyter/plotter_outer.py) and [Inner loop plotter](./jupyter/plotter_outer.py) functions

## Data
The format of the data is simple. Each file corresponds to a single run of the solver; each line corresponds to one outer loop (finding an augmenting step to apply according to the augmentation strategy), each number on a line corresponds to one run of the (AugILP) and is exactly the objective if this step would have been applied. The directory structure then corresponds to the various values of g1, augmentation strategy etc. (the directory names are self-explanatory).

- [Scheduling experimental data](./jupyter/sched/)
- [Closest string experimental data](./jupyter/csp/)

## Plots

We provide plots generated from data collected for scheduling instances with parameters m=15, lengths = {2,3,13,35}, weights = {6,13,2,1}, S=2000, L=10000, and slack ratio values 0.45, 0.50, 0.55, 0.60, 0.65.
We have tested the augmentation strategies "best step", "2-apx best step", "5-apx best step" and "any step".
The values of the parameter g1 were {2,3,4,5,6,7,8,9,10,12,14,16,18,20,23,26,29,33,40,45,50,60,75,90,110}.

### Outer loop plots
- [Ratio 0.45](./jupyter/graphing-longer_jobs-45-outer.ipynb)
- [Ratio 0.50](./jupyter/graphing-longer_jobs-50-outer.ipynb)
- [Ratio 0.55](./jupyter/graphing-longer_jobs-55-outer.ipynb)
- [Ratio 0.60](./jupyter/graphing-longer_jobs-60-outer.ipynb)
- [Ratio 0.65](./jupyter/graphing-longer_jobs-65-outer.ipynb)


### Inner loop plots
- [Ratio 0.45](./jupyter/graphing-longer_jobs-45-inner.ipynb)
- [Ratio 0.50](./jupyter/graphing-longer_jobs-50-inner.ipynb)
- [Ratio 0.55](./jupyter/graphing-longer_jobs-55-inner.ipynb)
- [Ratio 0.60](./jupyter/graphing-longer_jobs-60-inner.ipynb)
- [Ratio 0.65](./jupyter/graphing-longer_jobs-65-inner.ipynb)
