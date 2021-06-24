%%% Scope: Script for computing whole-cortex reconfiguration scores;
%%%         script contains filtering of connections and computes 
%%%         reconfiguration between connections that remain after the
%%%         filtering step
%%% Author: Jonas Thiele
%%% Date: 24.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% load parameters for analysis
load init_parameters

% load FCs (created in get_FCs)
load FCStatic_combined.mat

% load sample informations (created in get_samples)
load subjects_samples.mat

% load data from subjects (created in select_subjects)
behavioralData = readtable('data_beh_sel_subjects.csv');
subjects = behavioralData.Subject;

% load g-scores from subjects
load g_score_sel

% load confounds_sel for subjects
load confounds_sel.mat

% combinations of state pairs
combis_states = nchoosek((1:nStates),2);

% inizialise arrays for storing reconfiguration scores
dist_cosine_join = zeros(size(combis_states,1),length(subjects));
dist_corr_join = zeros(size(combis_states,1),length(subjects));
dist_L1_bin_join = zeros(size(combis_states,1),length(subjects));

for k=1:nSamples
    
    %% #1 Determine withheld sample and its parameters
    
    % all subjects but 1 sample (left out sample = kth sample)
    ind_subjects_filt = ~ismember(subjects, subjects_samples(k).IDs);
    subjects_filt = subjects(ind_subjects_filt);
    
    % FCs, g_scores and confounds of withheld sample
    FC_filt = FCStatic_combined(:,:,ind_subjects_filt,:);
    g_score_filt = g_score_sel(ind_subjects_filt);
    confounds_filt = confounds_sel(ind_subjects_filt,:);
    
    %% #2 Compute edges that are consitantly correlated to intelligence from withheld sample
    % filter for filtering connections in respect to their
    % consistency in their relation with intelligence
    % with intelligence
    
    % creating a mask for upper triangle of matrix
    mask_triu = triu(ones(nNodes));
    mask_triu = mask_triu - diag(diag(mask_triu));
    
    [nNodes, ~, nSubjects, ~] = size(FC_filt); 
    
    mask_pos = zeros(nNodes*nNodes,1);
    mask_neg = zeros(nNodes*nNodes,1);
    
    % collect all connections that are correlated with intelligence with 
    % p smaller threshold over all states
    for c = 1:nStates 
    
        FC_state = squeeze(FC_filt(:,:,:,c));
        FC_state = reshape(FC_state, nNodes*nNodes, nSubjects);
        [rho_fc, p_fc] = partialcorr(g_score_filt,FC_state',confounds_filt,'type','spearman');
        index_p = find(abs(p_fc) < p_thresh_whole_brain);

        mask_p = zeros(nNodes*nNodes,1);
        mask_p(index_p) = 1;

        mask_p_triu = mask_p.*mask_triu(:);
        
        % divide connections according to the directions of their correlations
        % with intelligence (positively or negatively) correlated
        mask_neg(((mask_p_triu ==1).* (rho_fc<0)')==1)=1;
        mask_pos(((mask_p_triu ==1).* (rho_fc>0)')==1)=1;
    
    end
    % mask of edges that are consistetly positively
    % correlated with intelligence across all scans
    mask_pos_fil = (mask_pos >= 1 & mask_neg == 0); 
    
    % mask of edges that are consistetly negatively
    % correlated with intelligence across all scans
    mask_neg_fil = (mask_pos == 0 & mask_neg >= 1); 
     
    %% #3 Compute filter mask from withheld sample
    
    for c = 1:length(combis_states)
        % states to be compared with 
        state1 = combis_states(c,1);
        state2 = combis_states(c,2);
        
        % compute edges positively correlating with intelligence 
        corr_mask_scan1_pos = get_correlation_mask(FC_filt, g_score_filt, confounds_filt, p_thresh_whole_brain, state1, 1);
        % filtering out inconsistant edges (mask computed above)
        mask_scan1_pos = corr_mask_scan1_pos.*mask_pos_fil;
        
        % compute edges negatively correlating with intelligence 
        corr_mask_scan1_neg = get_correlation_mask(FC_filt, g_score_filt, confounds_filt, p_thresh_whole_brain, state1, -1);
        mask_scan1_neg = corr_mask_scan1_neg.*mask_neg_fil;
    
        corr_mask_scan2_pos = get_correlation_mask(FC_filt, g_score_filt, confounds_filt, p_thresh_whole_brain, state2, 1);
        mask_scan2_pos = corr_mask_scan2_pos.*mask_pos_fil;
    
        corr_mask_scan2_neg = get_correlation_mask(FC_filt, g_score_filt, confounds_filt, p_thresh_whole_brain, state2, -1);
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
    
    mask_filter = mask_inter_all;
    % mask_filter = mask_triu_all; %%only to add if all edges should be used for calculating reconfiguration
    
    %% #4 Determine left out sample and its parameters
    
    subjects_ana = subjects_samples(k).IDs;
    ind_subjects_ana = ismember(subjects, subjects_ana);
    FC_ana = FCStatic_combined(:,:,ind_subjects_ana,:);
    g_score_ana = g_score_sel(ind_subjects_ana);
    confounds_ana = confounds_sel(ind_subjects_ana,:);
    
    %% #5 Compute reconfiguration
    [~, ~, nSubjects, ~] = size(FC_ana);
    
    dist_cosine = zeros(length(combis_states),nSubjects);
    dist_corr = zeros(length(combis_states),nSubjects);  
    dist_L1_bin = zeros(length(combis_states),nSubjects);

    for c = 1:length(combis_states)
    
        % states to be compared with
        state1 = combis_states(c,1);
        state2 = combis_states(c,2);

        for s=1:nSubjects % loop over subjects

            FCi = FC_ana(:,:,s,state1); 
            FC_filt_1 = FCi(:); % reshape FC for filtering
            % filter FC of scan with respective filter mask
            FC_filt_1 = FC_filt_1(mask_filter(:,c)==1); 

            FCi = FC_ana(:,:,s,state2); 
            FC_filt_2 = FCi(:);
            FC_filt_2 = FC_filt_2(mask_filter(:,c)==1);

            % calculate reconfiguration scores between filtered FCs of two
            % scans
            dist_cosine(c,s) = pdist([FC_filt_1';FC_filt_2'], 'cosine'); 
            dist_corr(c,s) = corr(FC_filt_1,FC_filt_2,'type','pearson'); 

            % compute bi-partitioned FCs and compare them
            FC_filt_1_bin = FC_filt_1 > 0;
            FC_filt_2_bin = FC_filt_2 > 0;

            dist_L1_bin(c,s) = norm(FC_filt_1_bin-FC_filt_2_bin,1) /sum(mask_filter(:,c));
        end
    end
    
    % add reconfiguration scores of the k-th sample to measures of all
    % subjects
    dist_cosine_join(:,ind_subjects_ana) = dist_cosine;
    dist_corr_join(:,ind_subjects_ana) = dist_corr;
    dist_L1_bin_join(:,ind_subjects_ana) = dist_L1_bin;
    
end

save('reconfiguration_scores.mat','dist_cosine_join','dist_corr_join','dist_L1_bin_join')
