# About these files

- .tar.bz2 files contain necessary experimental data (although paths in notebooks need to be changed to fit)
- [JEA.ipynb](./JEA.ipynb) is the main notebook containing most of necessary explanation
- [JEA_aux.ipynb](./JEA_aux.ipynb) is an auxiliary notebook for some data massaging
  (stripping pickles, cleaning full pickles for the purpose of publication etc.)
  
## How to generate experimental data and run the Jupyter notebook

- modify the path to the `sage` executable in [nfold_sched_tester.sage](./nfold_sched_tester.sage) to match
  wherever your installation of sage is
- install Gurobi and make it play nice with Sage
- create a virtualenv and install into it the following packages
  - `numpy`
  - `pandas`
  - `seaborn`
  - `jupyter`
  - (optionally) `jupyterlab` (simply because it's awesome.)
- activate the virtualenv

Then, use the [nfold_sched_tester.sage](./nfold_sched_tester.sage) script as documented
(e.g. just run it with `--help` and it will explain itself).

Specifically, the data sets demonstrated in the paper have been computed using the commands:

**Scheduling data set:** `./nfold_sched_tester.sage --logdir 31012019 --machines 10 20 30 40 50 60 --count_for_each_p 1 --slacks 0.6 0.7 --p_s 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 --number_job_types 4 5 6 --gammas log2 --gc 5 10 15 20 25 30 40 50 75 100 150 250 500 750 1000`

**Closest String data set:** `./nfold_sched_tester.sage.py --logdir logs2 --instance_type cs`

Beware that the computation of the Closest String data set takes a very long time because the dimensions of the generated instances are quite large.

## Additional explanation of plots

One might be wondering about the spikes in Figures 1 and 3.
Those represent the quality of augmentations found for various step lengths with respect to a fixed incumbent solution.
Specifically, we observe that typically, the best augmenting step is found neither for the shortest nor the longest step length, which accounts for the curve first having a "dip" in the middle, when viewed in the context of one outer loop iteration.

Some other spikes and fluctuations can be seen in Figures 8 and 9, which come from nondeterminism inherent in the algorithms and randomness in the data set.