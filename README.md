# Modeling the translation kinetics after mRNA transfection

This repository accompanies the article 

*Pieschner, Hasenauer, Fuchs (2021). "Identifiability analysis for models of the translation kinetics after mRNA transfection"*.

**Abstract**

Mechanistic models are a powerful tool to gain insights into biological processes. The parameters of such models, e.g. kinetic rate constants, usually cannot be measured directly but need to be inferred from experimental data. In this article, we study dynamical models of the translation kinetics after mRNA transfection and analyze their parameter identifiability. Previous studies have considered  ordinary differential equation (ODE) models of the process, and here we formulate a stochastic differential equation (SDE) model. For both model types, we consider structural identifiability based on the model equations and practical identifiability based on simulated as well as experimental data and find that the SDE model provides better parameter identifiability than the ODE model. Moreover, our analysis shows that even for those parameters of the ODE model that are considered to be identifiable, the obtained estimates are sometimes unreliable. Overall, our study clearly demonstrates the relevance of considering different modeling approaches and that stochastic models can provide more reliable and informative results.

In this project, we use simulated data generated as described below and experimental data that has been published in Fröhlich et al (2018, 
https://doi.org/10.1038/s41540-018-0079-7, see `data/experimental_data/README.md` for further details).

### Repository content
This repository contains all relevant code for the data analysis.

The computational environment is described in the Docker file in `container_image/dockerfile_r3.6.2_rstan_rmd`.

#### Code related to structural identifiability analysis
The model and output files for the structural identifiability analysis with DAISY are contained in the folder `Identifiability_analysis_with_DAISY`. See `Identifiability_analysis_with_DAISY/README.md` for further details.

The figures for the simulation based assessment of parameter influence were generated in `R` with `R_code/figures_and_tables/figs_simulate_trajectories_to_study_identifiability_with_param_from_Froehlich.R`.

#### Code related to practical identifiability analysis

We use the open source software `Stan` to sample from the posterior distributions of the ODE or SDE model given the experimental or simulated data, respectively. The Stan models are defined in the `Stan_model_code` folder.

We use `Stan` through its interface `rstan` in `R`. All `R` code for generating simulated data, sampling from the posterior distribution, and aggregating the results is contained in the `R_code` folder.  
Moreover, the `R` code to generate figures and table for the article is contained in `R_code/figures_and_tables`.

Most of the calculations were performed on a HPC cluster managed by SLURM. The bash scripts to submit the jobs are contained in the folder `bash_scripts`.

Additional summaries of all sampling results and for some individual trajectories are contained in the folder `R_markdown`.
