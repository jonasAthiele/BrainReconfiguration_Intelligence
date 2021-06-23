
%%% Scope:  Script for defining parameters for analysis
%%% Author: Jonas Thiele
%%% Date:   16.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% folder and filename of .nii data (fMRI data)
folderRest = 'regress_36pNS_schaefer200';
folderTask = 'taskregress_36pNS_schaefer200-yeo17';
filename_nii = 'schaefer200-yeo17.ptseries.nii';

% number of regions (nodes) of parcellation used, range [1 to nNodes]
% if range [>1 to nNodes] should be used, change code in get_FCs
nNodes = 200; % in total we have 255 nodes, 1:200 = cortical, 201:255 = subcortical and cerebellum 

% list of scans to be considered
scans = {'rfMRI_REST1_RL','rfMRI_REST1_LR','rfMRI_REST2_RL','rfMRI_REST2_LR','tfMRI_WM_RL','tfMRI_WM_LR',...
   'tfMRI_GAMBLING_RL', 'tfMRI_GAMBLING_LR', 'tfMRI_MOTOR_RL', 'tfMRI_MOTOR_LR', 'tfMRI_LANGUAGE_RL',...
   'tfMRI_LANGUAGE_LR', 'tfMRI_SOCIAL_RL', 'tfMRI_SOCIAL_LR', 'tfMRI_RELATIONAL_RL', 'tfMRI_RELATIONAL_LR',...
   'tfMRI_EMOTION_RL', 'tFMRI_EMOTION_LR'};

% indexes of resting state in scans
inds_rest = 1:4;

% mask for joining scans (done in get_FCs code)
mask_states = [1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8]; % which tasks belong together 4x rest, 2 per task.
nStates = max(mask_states); % number of scan conditions (fMRI states) after joining


% k-folds for cross-validated approach - k is equal to nSamples
nSamples = 10;

% threshold p-value for filtering connections of FCs
% FCs are correlated with g-score, only edges with p smaller p_tresh are
% kept
p_thresh_whole_brain = 0.1;
p_thresh_yeo7 = 0.1;
p_thresh_yeo17 = 0.2;

% threshold p-value for filtering connections of replication samples FCs
p_thresh_ref_mask = 0.01;

% alpha value for FDR correction
alpha_FDR = 0.05;

% number of functional networks for network-wise reconfiguration must be 7
% or 17
nNetworks = 7; % 7 or 17 ONLY!!

% save variables
save init_parameters

