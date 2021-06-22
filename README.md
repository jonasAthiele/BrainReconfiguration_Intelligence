# BrainReconfiguration_Intelligence

## 1. Scope
The repository contains scripts for the analysis used in the paper **"Multi-Task Brain Network Reconfiguration and its Inverse Association with General Intelligence"** coauthored Jonas A. Thiele, Joshua Faskowitz, Olaf Sporns, and Kirsten Hilger (doi: coming soon after publication). Herein, the relation between general intelligence and the reconfiguration of functional brain connectivity is analyzed.
The scripts in this repository can be used to replicate the analyses of the paper or more generally, to study associations between individual differences (e.g., in intelligence) and the reconfiguration of functional brain connectivity.
If you have questions or trouble with the scripts, feel free to contact me: jonas.thiele@uni-wuerzburg.de
## 2. Data
We used data provided by the Human Connectome project (1), funded by the National Institute of Health for providing data for our main sample analysis. Data from The Amsterdam Open MRI Collection (2) was used for the replication analyses (PIOP1 and PIOP2 sample).
All data used in the current study can be accessed online under: https://www.humanconnectome.org/study/hcp-young-adult (HCP), https://doi.org/10.18112/openneuro.ds002785.v2.0.0 (AOMIC-PIOP1), and https://doi.org/10.18112/openneuro.ds002790.v2.0.0 (AOMIC-PIOP2).
## 3. Preprocessing
We used the minimally preprocessed HCP fMRI data (1) and implemented further preprocessing comprising a nuisance regression strategy with 24 head motion parameters, eight mean signals from white matter and cerebrospinal fluid, and four global signals (3). For task data, basis-set task regressors (4) were used simultaneously with the nuisance regressors to remove mean task-evoked activations.
Code for the further preprocessing steps is available here: https://github.com/faskowit/app-fmri-2-mat
For the replication the data of The Amsterdam Open MRI Collection was downloaded in the minimal preprocessed form (using fMRIPrep version 1.4.1) and all further preprocessing followed the same regression steps as specified for the main sample.  
For all data, timeseries of neural activation were extracted from 200 nodes covering the entire cortex (5).
## 4. Structure and Script description
### Main analysis
For the analysis done in the paper the scripts should be run in the following order:
1.	`get_init_parameters` (Matlab)
  - Script for defining parameters for analysis --> can be adapted for modifying the analysis
  
  
2.	`select_subjects` (Matlab)
  - Script for the selection of subjects according to data availability and motion exclusion criteria
  
3.	`get_gFactor` (R)
  - Derive a g-factor from the performance scores of cognitive tests
  
4.	`get_samples` (Matlab)
  - Create samples with absence of family relations between subjects and an equal distribution of intelligence scores via stratified folds
  
5.	`get_FCs` (Matlab)
  - Reading preprocessed BOLD-signals, calculating FCs, joining FCs
  
6.	`get_reconfiguration_wholebrain` (Matlab)
  - Script for computing whole-brain reconfiguration scores, contains filtering of connections and computes reconfiguration between connections that remain after the filtering step
7.	`get_relation_recon_intell_wholebrain`
  - Script for computing correlations between reconfiguration scores and intelligence scores on a whole-cortex level


**and/or**


6.	`get_reconfiguration_networks` (Matlab)
  - Script for computing within and between network reconfiguration, contains filtering of connections and computes reconfiguration between connections that remain after the filtering step
7.	`get_relation_recon_intell_networks` (Matlab)
  - Script for computing correlations between reconfiguration scores and intelligence scores on a network level
  
8.	 `get_hcp_reference_mask` (Matlab)
  - Script for computing reference masks for edge-selection in the replication samples
  
### For replication analysis

For the replication analysis done in the paper the scripts should be run in the following order:

1.	`get_init_parameters_replication` (Matlab)
  - Script for defining parameters for analysis
2.	`select_subjects_replication` (Matlab)
  - Script for the selection of subjects according to data availability and motion exclusion criteria
3.	`get_FCs_replication` (Matlab)
  - Reading Preprocessed BOLD-signals, calculating FCs, joining FCs
4.	`get_reconfiguration_wholebrain_replication` (Matlab)
  - Script for computing whole-brain reconfiguration scores, contains filtering of conenctions and computes reconfiguration between connections that remain after filtering step
5.	`get_relation_recon_intell_wholebrain` (same as in main analysis)
  - Script for computing correlations between reconfiguration scores and intelligence scores on a whole-cortex level 


**and/or**


4.	`get_reconfiguration_networks_replication` (Matlab)
  - Script for computing within and between network reconfiguration, contains filtering of connections and computes reconfiguration between connections that remain after filtering step
5.	`get_relation_recon_intell_networks` (same as in main analysis)
  - Script for computing correlations between reconfiguration scores and intelligence scores on a network level

### Functions 

Functions used in the scripts can be found in the “Functions” folder. Some of the functions can be found elsewhere but are included here for convenience. Comments on the authorship of all the functions are provided within the functions.

## 5. Software requirements
-	Matlab 2021a
-	R version 4.0.2

## References
1.	D. C. Van Essen, et al., The WU-Minn Human Connectome Project: An overview. Neuroimage 80, 62–79 (2013).
2.	L. Snoek, et al., The Amsterdam Open MRI Collection, a set of multimodal MRI datasets for individual difference analyses. Sci. Data 8, 85 (2021).
3.	L. Parkes, B. Fulcher, M. Yücel, A. Fornito, An evaluation of the efficacy, reliability, and sensitivity of motion correction strategies for resting-state functional MRI. Neuroimage 171, 415–436 (2018).
4.	M. W. Cole, et al., Task activations produce spurious but systematic inflation of task functional connectivity estimates. Neuroimage 189, 1–18 (2019).
5.	A. Schaefer, et al., Local-Global Parcellation of the Human Cerebral Cortex from Intrinsic Functional Connectivity MRI. Cereb. Cortex 28, 3095–3114 (2018).

## Copyright
Copyright (cc) 2021 by Jonas Thiele


<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">Files of BrainReconfiguration_Intelligence</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="https://github.com/jonasAthiele/BrainReconfiguration_Intelligence" property="cc:attributionName" rel="cc:attributionURL">Jonas A. Thiele</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License</a>.
