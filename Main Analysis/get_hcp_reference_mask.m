
%%% Scope:  Script for computing reference masks for edge-selection 
%%%         in the replication samples
%%% Author: Jonas Thiele
%%% Date:   24.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% load parameters for analysis
load init_parameters

% load data from subjects (created in select_subjects)
behavioralData = readtable('data_beh_sel_subjects.csv');
subjects = behavioralData.Subject;

% load FCs (created in get_FCs)
load FCStatic_combined.mat

% load g-scoresfrom subjects
g_score_1186 = readtable('gfac_bi_explore_1186Subjects.csv').x;
subjects_1186 = readtable('subject_IDs_g_score_1186Subjects.csv').Subject;
ind_sel = ismember(subjects_1186, subjects);
g_score_sel = g_score_1186(ind_sel);

% load confounds_sel for subjects
load confounds_sel.mat


[nNodes, ~, nSubjects, nStates] = size(FCStatic_combined);

% mask for upper triangle of matrix without diagonal
mask_triu = triu(ones(nNodes));
mask_triu = mask_triu - diag(diag(mask_triu));
combis_states = nchoosek((1:nStates),2);

mask_neg = zeros(nNodes*nNodes,1);
mask_pos = zeros(nNodes*nNodes,1);

for c = 1:length(combis_states)
    
    % states to be compared with
    state1 = combis_states(c,1);
    state2 = combis_states(c,2);
   
    % relation between edge strenghts and intelligence
    FC_state1 = squeeze(FCStatic_combined(:,:,:,state1));
    FC_state1 = reshape(FC_state1, nNodes*nNodes, nSubjects);
    [rho_fc, p_fc] = partialcorr(g_score_sel,FC_state1',confounds_sel,'type','spearman');
    
    % find p-values of relations smaller p-threshold
    index_p = find(abs(p_fc) < p_thresh_ref_mask);
    
    % compute mask with p-values smaller p-threshold
    mask_p = zeros(nNodes*nNodes,1);
    mask_p(index_p) = 1;
    
    % multiply p-mask with mask of upper triangle - only use p-values in
    % upper triangle of matrix 
    mask_p_triu = mask_p.*mask_triu(:);
    
    % find correlations that are positive (mask_pos) and negative (mask_neg)
    mask_neg_1 = (mask_p_triu ==1).* (rho_fc<0)';
    mask_pos_1 = (mask_p_triu ==1).* (rho_fc>0)';
    
    % same as above for FCs of second state
    FC_state2 = squeeze(FCStatic_combined(:,:,:,state2));
    FC_state2 = reshape(FC_state2, nNodes*nNodes, nSubjects);
    [rho_fc, p_fc] = partialcorr(g_score_sel,FC_state2',confounds_sel,'type','spearman');
    index_p = find(abs(p_fc) < p_thresh_ref_mask);

    mask_p = zeros(nNodes*nNodes,1);
    mask_p(index_p) = 1;

    mask_p_triu = mask_p.*mask_triu(:);

    mask_neg_2 = (mask_p_triu ==1).* (rho_fc<0)';
    mask_pos_2 = (mask_p_triu ==1).* (rho_fc>0)';
    
    % compute intersecting masks of state 1 and 2
    mask_neg_com = mask_neg_1.*mask_neg_2;
    mask_pos_com = mask_pos_1.*mask_pos_2;

    % add edges to mask of all state comparisons 
    mask_neg(mask_neg_com==1)=1;
    mask_pos(mask_pos_com==1)=1;

end

% mask of edges that are consistetly positively
% correlated with intelligence across all state comparisons
mask_pos_fil = (mask_pos >= 1 & mask_neg == 0); 

% mask of edges that are consistetly negatively
% correlated with intelligence across all state comparisons
mask_neg_fil = (mask_pos == 0 & mask_neg >= 1); 
    
save('reference_mask_hcp', 'mask_pos_fil', 'mask_neg_fil');    