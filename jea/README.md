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
