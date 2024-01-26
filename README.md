# CL-DeePC and its equivalence with CL-SPC
This repository contains the necessary components of the Automatica article titled "Closed-loop Data-Enabled Predictive Control and its equivalence with Closed-loop Subspace Predictive Control".

## Project organization
- PG = project-generated
- HW = human-writable
- RO = read only
```
.
├── .gitignore
├── CITATION.md
├── LICENSE.md
├── README.md
├── bin                         <- Compiled and external code, ignored by git (PG)
│   ├── external/Casadi         <- When running simulations on a cluster, place compiled CasADi code with IPOPT solver here. (RO)
|   └── model_Favoreel1999.m    
├── data                        <- All project data, ignored by git
│   ├── processed               <- The final, canonical data sets for modeling. (PG)
│   ├── raw                     <- The original, immutable data dump. (RO)
│   └── temp                    <- Intermediate data that has been transformed. (PG)
├── docs                        <- Documentation notebook for users (HW)
│   └── manuscript              <- Overleaf code (HW)
├── results         
│   └── figures                 <- Figures for the manuscript or reports (PG)
└── src                         <- Source code for this project (HW)

```
## Configuration & use:
The generated data was obtained by varying either one of three parameters:
- Nbar: number of samples used
- p=f:  window sizes
- Re:   innovation noise variance
The functions to generate the relevant data and a sample script for submission of multiple jobs to a Slurm cluster can be found in the corresponding directories src/d_{Nbar,pf,Re}.

The data that is used in the article can be obtained in one of two ways:
1. Download it from [link to Zenodo] and place it in the data\raw directory.
2. Generate it data. To do this follow the below steps:
    Step 1: Generate the the necessary data using the Slurm submission script (or make a similar Matlab script) to run different simulation cases with the functions 'varying_{...}.m' found in the relevant aformentioned directory.
    Step 2: The generated data should automatically be placed in a descriptively named folder in the data/raw directory. Go to this new folder with the new data and run 'src/dirs2res.m'.

To generate the figures in the paper run the script 'src/make_figures.m'. The figures will be placed inside the results/figures directory.

## License

This project is licensed under the terms of the [MIT License](/LICENSE.md)
