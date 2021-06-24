
%%% Scope:  Script for computing correlations between reconfiguration 
%%%         scores and intelligence scores on a whole-cortex level
%%% Author: Jonas Thiele
%%% Date:   26.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% load parameters for analysis
load init_parameters

% load data from subjects (created in select_subjects)
behavioralData = readtable('data_beh_sel_subjects.csv');
subjects = behavioralData.Subject;

% load g-scores from subjects
load g_score_sel

% load confounds_sel for subjects
load confounds_sel.mat

% load reconfiguration scores
load reconfiguration_scores

%% Compute relation between total rest-task and task-task reconfiguration and intelligence

% distances between rest and task
[rho_mean_dist_rest_task, p_mean_dist_rest_task] = partialcorr(g_score_sel,mean(dist_cosine_join(1:nStates-1,:))',confounds_sel,'type','spearman');
[rho_mean_corr_rest_task, p_mean_corr_rest_task] = partialcorr(g_score_sel,mean(dist_corr_join(1:nStates-1,:))',confounds_sel,'type','spearman');
[rho_mean_L1bin_rest_task, p_mean_L1bin_rest_task] = partialcorr(g_score_sel,mean(dist_L1_bin_join(1:nStates-1,:))',confounds_sel,'type','spearman');

% regress out confounds fromn reconfiguration scores for scatter plot
y = mean(dist_cosine_join(1:nStates-1,:))';
X = [ones(size(y)), confounds_sel];
[~,~,residuals_rest_task] = regress(y,X);

g_score_vis = g_score_sel;
res_dist_cosine_rest_task = normalize(residuals_rest_task); % normalized residuals

%%%% subject 243 excluded for visualization
% g_score_vis(243) = []; 
% res_dist_cosine_rest_task(243)=[];

figure()
subplot(1,2,1)
scat1=scatter(g_score_vis,res_dist_cosine_rest_task,'k.');
scat1.SizeData = 65;
h=lsline;
h.LineWidth = 2;
h.Color = 'k';
xlabel('g-score') 
ylabel('reconfiguration std. res.') 
title('rest-task reconfiguration')

% distances between tasks
[rho_mean_dist_task_task, p_mean_dist_task_task] = partialcorr(g_score_sel,mean(dist_cosine_join(nStates:end,:))',confounds_sel,'type','spearman');
[rho_mean_corr_task_task, p_mean_corr_task_task] = partialcorr(g_score_sel,mean(dist_corr_join(nStates:end,:))',confounds_sel,'type','spearman');
[rho_mean_L1bin_task_task, p_mean_L1bin_task_task] = partialcorr(g_score_sel,mean(dist_L1_bin_join(nStates:end,:))',confounds_sel,'type','spearman');

% regress out confounds fromn reconfiguration scores for scatter plot
y = mean(dist_cosine_join(nStates:end,:))';
X = [ones(size(y)), confounds_sel];
[~,~,residuals_task_task] = regress(y,X);

res_dist_cosine_task_task = normalize(residuals_task_task); 

%%%% subject excluded from plot, seems to be an outlier
% res_dist_cosine_task_task(243)=[]; 

subplot(1,2,2)
scat2 = scatter(g_score_vis,res_dist_cosine_task_task,'k.');
scat2.SizeData = 65;
h = lsline;
h.LineWidth = 2;
h.Color = 'k';
xlabel('g-score') 
ylabel('reconfiguration std. res.') 
title('task-task reconfiguration')

%% Compute relation between reconfiguration and intelligence for each pair of states

range_plot = [-0.23, 0.23]; % range of values shown in colorbar of plots

% cosine distance
[rho_cosdist_pairs, p_cosdist_pairs] = partialcorr(g_score_sel,dist_cosine_join',confounds_sel,'type','spearman');
figure()
subplot(1,3,1)
plot_rho_p_fdrCorrected(rho_cosdist_pairs, p_cosdist_pairs, alpha_FDR, nStates, range_plot, 'cos dist - intelligence')

% FDR adjusted p-values
[~, ~, ~, p_cosdist_pairs_adj] = fdr_bh(p_cosdist_pairs(isfinite(p_cosdist_pairs)), alpha_FDR);

% correlation
[rho_corr_pairs, p_corr_pairs] = partialcorr(g_score_sel, dist_corr_join',confounds_sel,'type','spearman');
subplot(1,3,2)
plot_rho_p_fdrCorrected(rho_corr_pairs, p_corr_pairs, alpha_FDR, nStates, range_plot, 'corr - intelligence')

% bi-partition Manhatten distance
[rho_L1_bin_pairs, p_L1_bin_pairs] = partialcorr(g_score_sel, dist_L1_bin_join',confounds_sel,'type','spearman');
subplot(1,3,3)
plot_rho_p_fdrCorrected(rho_L1_bin_pairs, p_L1_bin_pairs, alpha_FDR, nStates, range_plot, 'L1 bin dist - intelligence')

%% Relation between total reconfiguration for each state and intelligence

combis_states = nchoosek(1:nStates,2);

% compute a total (average) reconfiguration score for each state
% average over state comparisons the respective state is involved in
for c=1:nStates
    
    ind = combis_states(:,1) == c | combis_states(:,2) == c;
    if c>1
        ind(1:nStates-1)=0;
    end
    
    distcos_total_states(:,c) = mean(dist_cosine_join(ind,:));
    [rho_total_states(c), p_total_states(c)] = partialcorr(g_score_sel, mean(dist_cosine_join(ind,:))',confounds_sel,'type','spearman');
    
end


figure()
bar(rho_total_states,'FaceColor',[0.7,0.7,0.7])
ylim([-0.3, -0.1])
xlabel('states')
ylabel('rho')
title('rho: state specific reconfiguration - intelligence')

% check which differenes in the relations (total reconfiguration -
% intelligence) between the states are significant
for c = 1:length(combis_states)
    
    t1 = combis_states(c,1);
    t2 = combis_states(c,2);
    dist_t1 = distcos_total_states(:,t1);
    dist_t2 = distcos_total_states(:,t2);
       
    
    X = [ones(size(y)), confounds_sel];
    
    [~,~,residuals_t1] = regress(dist_t1,X);
    [~,~,residuals_t2] = regress(dist_t2,X);
    
    [p_t(c), T2(c), df(c), r_jk(c), r_jh(c), r_kh(c)] = ...
        r_test_paired(g_score_sel, normalize(residuals_t1)',normalize(residuals_t2)',0);
 
end

M_tasks = NaN(nStates);

for c=1:length(combis_states)
    
    el1 = combis_states(c,1);
    el2 = combis_states(c,2);
    M_tasks(el1,el2) = p_t(c);
   
end

figure()
h=imagesc(M_tasks);
colorbar
line(repmat([1.5:1:nStates-0.5],2,1),repmat([0;nStates+1], 1, nStates-1), 'Color', 'black'); % vertical
line(repmat([0;nStates+1], 1, nStates-1), repmat([1.5:1:nStates-0.5],2,1), 'Color', 'black'); % horizontal
set(h, 'AlphaData', ~isnan(M_tasks));  
xlabel('states')
ylabel('states')
title('p-value of differences in corr (g - reconfiguration) between states')


