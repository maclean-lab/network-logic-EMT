# Transition paths across the EMT landscape are dictated by network logic

### Citation
Please cite the [paper](https://www.biorxiv.org/content/10.1101/2024.12.03.626660v1): 

A Dey, AL MacLean (2024), Transition paths across the EMT landscape are dictated by network logic. bioRxiv, 10.1101/2024.12.03.626660. 

### Overview 
This repository contains MATLAB code for analyzing the effects of network logic (AND vs OR) on multistable models of EMT. We analyze tristability via estimation of a potential energy function that permits calculation of the bifurcation landscapes, and analyze model behaviors via statistical analyses and  stochastic simulations. 

Requirements
- MATLAB (tested on version R2024a)  
- All parameter values are included within respective code files unless otherwise specified
- Generated parameter files are required for certain analyses as detailed below.

### Project Contents and Usage 
The following describes the contents of this repository, which can be used to reproduce all of the figures of the paper. The code required to reproduce each figure is contained within the folder of that name. To recapitulate the results from the paper. 

1. Run the parameter search code first (`CodeMain_RandParSearch_AND.m` and `CodeMain_RandParSearch_OR.m`).
2. Run perturbation analysis code (`CodeMain_PertAND.m` and `CodeMain_PertOR.m`).
3. Run stochastic simulations.
4. Generate figures using the respective figure codes with required parameter files.


### Figure 1
- `Code_Fig1B.m` - Generates the multistable landscapes plotted in Fig. 1B
  - Self-contained code with embedded parameters

### Figure 2 (AND Model)
- `Code_Fig2AC.m` - Generates Figures 2A and 2C
  - Modify SA0 (-25%, -50%) for Fig 2C
  - Modify SB0 (-25%, -50%) for Fig 2A
- `Code_Fig2B.m` - Generates Figure 2B
- `Code_Fig2D.m` - Generates Figure 2D
- `CodeMain_RandParSearch_AND.m` - Random parameter space search for tristability
  - Outputs stored in `param_SameSA0SB0_AND.mat`
- `CodeMain_PertAND.m` - Perturbation analysis of AND models
  - Requires `param_SameSA0SB0_AND.mat`
  - Outputs stored in `param_ThrePert_SA0SB0_25pc_AND.mat`
- `Code_Fig2EF.m` - Generates Figures 2E and 2F
  - Requires both parameter files above

### Figure 3 (OR Model)
- `Code_Fig3AC.m` - Generates Figures 3A and 3C
  - Modify SA0 (-25%, -50%) for Fig 3C
  - Modify SB0 (-25%, -50%) for Fig 3A
- `Code_Fig3B.m` - Generates Figure 3B
- `Code_Fig3D.m` - Generates Figure 3D
- `CodeMain_RandParSearch_OR.m` - Random parameter space search for tristability
  - Outputs stored in `param_SameSA0SB0_OR.mat`
- `CodeMain_PertOR.m` - Perturbation analysis of OR models
  - Requires `param_SameSA0SB0_OR.mat`
  - Outputs stored in `param_ThrePert_SA0SB0_25pc_OR.mat`
- `Code_Fig3EF.m` - Generates Figures 3E and 3F
  - Requires both parameter files above

### Figure 4
- `Fig4_BD_OR_25pcPert.m` - OR model analysis and generates Figs. 4B and 4D
  - Requires `param_SameSA0SB0_OR.mat` and `param_ThrePert_SA0SB0_25pc_OR.mat`
- `Fig4_AC_AND_25pcPert.m` - AND model analysis and generates Figs. 4A and 4C
  - Requires `param_SameSA0SB0_AND.mat` and `param_ThrePert_SA0SB0_25pc_AND.mat`
- `Fig4_E_AND.m` - Generates Fig. 4E for AND model
  - Requires `param_SameSA0SB0_AND.mat` and `param_ThrePert_SA0SB0_25pc_AND.mat`
- `Fig4_F_OR.m` - Generates Fig. 4F for OR model
  - Requires `param_SameSA0SB0_OR.mat` and `param_ThrePert_SA0SB0_25pc_OR.mat`

### Figure 5 (Stochastic model simulation)
- `AND_Multi_EulerMaryama.m` - SDE simulations for AND model
  - Outputs:
    - Unperturbed: `stochEM_sigma0o01_V1_Multi_Unpert_AND.mat`
    - 25% indirect: `stochEM_sigma0o01_V1_Multi_25pcInd_AND.mat`
    - 25% direct: `stochEM_sigma0o01_V1_Multi_25pcDir_AND.mat`
- `OR_Multi_EulerMaryama.m` - SDE simulations for OR model
  - Outputs:
    - Unperturbed: `stochEM_sigma0o01_V1_Multi_Unpert_OR.mat`
    - 25% indirect: `stochEM_sigma0o01_V1_Multi_25pcInd_OR.mat`
    - 25% direct: `stochEM_sigma0o01_V1_Multi_25pcDir_OR.mat`
- `Code_Fig5A.m` - Generates Figure 5A
  - Requires stochastic simulation outputs
- `Code_Fig5BCD.m` - Generates Figures 5B, 5C, and 5D
  - Requires stochastic simulation outputs

