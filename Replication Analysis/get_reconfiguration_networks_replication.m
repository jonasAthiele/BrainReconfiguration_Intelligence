
%%% Scope:  Script for computing within and between network 
%%%         reconfiguration; script contains filtering of connections 
%%%         and computes reconfiguration between connections that remain
%%%         after the filtering step
%%% Author: Jonas Thiele
%%% Date:   24.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% initialize parameters for analysis
load init_parameters

% load FCs (created in get_FCs)
load FCStatic_combined.mat


% load data from subjects (created in select_subjects)
behavioralData = readtable('data_beh_sel_subjects.csv');
subjects = behavioralData.Subject;

% load g-scores from subjects
load g_score_sel.mat

% load confounds_sel for subjects
load confounds_sel.mat

% combinations of scan pairs
combis_states = nchoosek((1:nStates),2);

% load assignments of nodes to networks 
if nNetworks == 7
    load Schaf200Yeo7_PIOP
    node_assignment = Schaf200Yeo7_PIOP;
elseif nNetworks == 17
    load Schaf200Yeo17_PIOP
    node_assignment = Schaf200Yeo17_PIOP;
else
    error('Wrong number of networks: nNetwork must be either 7 or 17, check init_parameters')
end

% create masks for the connections of within and between network
% combinations
mask_networks = get_mask_networks(node_assignment);
mask_networks = reshape(mask_networks,nNodes*nNodes,[]);

    
%% #1 Determine sample and its parameters
% renaming for compliance with scripts for HCP data
FC_ana = FCStatic_combined; 
g_score_ana = g_score_sel;
confounds_ana = confounds_sel;
    
%% #2 Compute filter mask 

% load reference mask derived from HCP sample 
% mask_pos_fil (edges correlating positively with intelligence)
% mask_neg_fil (edges correlating negatively with intelligence)
load reference_mask_hcp.mat

% creating a mask for upper triangle of the matrix
mask_triu = triu(ones(nNodes));
mask_triu = mask_triu - diag(diag(mask_triu));

for c = 1:length(combis_states)
    
    % states to be comared with 
    state1 = combis_states(c,1);
    state2 = combis_states(c,2);

    % compute edges positively correlating with intelligence 
    corr_mask_scan1_pos = get_correlation_mask(FC_ana, g_score_ana, confounds_ana, p_thresh, state1, 1);
    % filtering out inconsistant edges (mask computed above)
    mask_scan1_pos = corr_mask_scan1_pos.*mask_pos_fil;

    % compute edges negatively correlating with intelligence 
    corr_mask_scan1_neg = get_correlation_mask(FC_ana, g_score_ana, confounds_ana, p_thresh, state1, -1);
    mask_scan1_neg = corr_mask_scan1_neg.*mask_neg_fil;
    
    corr_mask_scan2_pos = get_correlation_mask(FC_ana, g_score_ana, confounds_ana, p_thresh, state2, 1);
    mask_scan2_pos = corr_mask_scan2_pos.*mask_pos_fil;

    corr_mask_scan2_neg = get_correlation_mask(FC_ana, g_score_ana, confounds_ana, p_thresh, state2, -1);
    mask_scan2_neg = corr_mask_scan2_neg.*mask_neg_fil;

    % combining connections correlating positively or negatively with
    % intelligence
    mask_scan1_combi = (mask_scan1_pos + mask_scan1_neg);
    mask_scan2_combi = (mask_scan2_pos + mask_scan2_neg);
    
    % interception between masks of both scans to compare with
    mask_inter = mask_scan1_combi.*mask_scan2_combi.*mask_triu(:);
    mask_inter(mask_inter>0) = 1;

    % store masks for each scan combination
    mask_inter_all(:,c) = mask_inter;
    mask_triu_all(:,c) = mask_triu(:);

end

% mask_filter = mask_triu_all; %%only to add if all edges should be used for calculating reconfiguration
mask_filter = mask_inter_all;



%% #3 Compute reconfiguration
[~, ~, nSubjects, ~] = size(FC_ana);

dist_cosine = zeros(size(mask_networks,2), length(combis_states), nSubjects);
dist_corr = zeros(size(mask_networks,2), length(combis_states), nSubjects);  
dist_L1_bin = zeros(size(mask_networks,2), length(combis_states), nSubjects);
for n = 1:size(mask_networks,2)
    for c = 1:length(combis_states)
        
        % states to be compared with
        state1 = combis_states(c,1);
        state2 = combis_states(c,2);

        % intersection between filter and network masks
        mask_filter_yeo = mask_networks(:,n).*mask_filter(:,c);

        for s=1:nSubjects % loop over subjects

        FCi = FC_ana(:,:,s,state1); 
        FC_filt_1 = FCi(:); % reshape FC for filtering
        % filter FC of scan with respective filter mask
        FC_filt_1 = FC_filt_1(mask_filter_yeo==1); 

        FCi = FC_ana(:,:,s,state2); 
        FC_filt_2 = FCi(:);
        FC_filt_2 = FC_filt_2(mask_filter_yeo==1);
        
        % check if no connection was kept after filtering
        % if so, set respective value NaN
        if sum(mask_filter_yeo)==0

            dist_cosine(n,c,s)= NaN;
            dist_corr(n,c,s) = NaN;
            dist_L1_bin(n,c,s) = NaN;


        else

            % calculate reconfiguration scores between filtered FCs of two
            % scans
            dist_cosine(n,c,s) = pdist([FC_filt_1';FC_filt_2'], 'cosine'); 
            dist_corr(n,c,s) = corr(FC_filt_1,FC_filt_2,'type','pearson'); 

            % compute bi-partitioned FCs and compare them
            FC_filt_1_bin = FC_filt_1 > 0;
            FC_filt_2_bin = FC_filt_2 > 0;

            dist_L1_bin(n,c,s) = norm(FC_filt_1_bin-FC_filt_2_bin,1)/sum(mask_filter(:,c));
        end

        end
    end
end

% renaming for compliance with scripts for HCP data
dist_cosine_join = dist_cosine; 
dist_corr_join = dist_corr;
dist_L1_bin_join = dist_L1_bin;
    
save('reconfiguration_scores_yeo.mat','dist_cosine_join','dist_corr_join','dist_L1_bin_join')
