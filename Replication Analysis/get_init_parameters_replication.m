
%%% Scope:  Script for defining parameters for replication analysis 
%%% Author: Jonas Thiele
%%% Date:   16.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% filename of (fMRI data)
filename_fMRI = 'out_schaefer200-yeo17_timeseries.hdf5';

% number of regions (nodes) of parcellation used 1 to nNodes
% if range >1 to nNodes should be used, change code in get_FCs
nNodes = 200; % in total we have 216 nodes, 1:202 = cortical with two background regions 203:216 subcortical nodes
background_regions = [1,102];
% background regions for Schaefer 100: [1,52] 

% list of scans to be considered
%%% scans PIOP1
scans = {'task-restingstate_acq-mb3','task-workingmemory_acq-seq',...
    'task-anticipation_acq-seq','task-emomatching_acq-seq',...
    'task-faces_acq-mb3','task-gstroop_acq-seq'};

%%% scans PIOP2
% scans = {'task-restingstate_acq-seq','task-workingmemory_acq-seq',...
%     'task-emomatching_acq-seq','task-stopsignal_acq-seq'};


nStates = length(scans); % number of scan conditions (fMRI states) after joining

% alpha value for FDR correction
alpha_FDR = 0.05;

% number of functional networks for network-wise reconfiguration must be 7
% or 17
nNetworks = 7; % 7 or 17 ONLY!!

% p-threshold for filtering connection, here set to 1 --> all connections
% (edges) are selected, connections are filtered with the
% reference_mask_hcp (filter mask derived from the HCP-sample)
p_thresh = 1;

save init_parameters