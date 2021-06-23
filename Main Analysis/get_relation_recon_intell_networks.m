
%%% Scope:  Script for computing correlations between reconfiguration 
%%%         scores and intelligence scores on a network level
%%% Author: Jonas Thiele
%%% Date:   16.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all

% load parameters for analysis
load init_parameters

% load data from subjects (created in select_subjects)
behavioralData = readtable('data_beh_sel_subjects.csv');
subjects = behavioralData.Subject;

% load g-scoresfrom subjects
load g_score_sel

% load confounds for subjects
load confounds_sel.mat

% load reconfiguration scores
load reconfiguration_scores_yeo

% reference mask with all within and between network combinations
[ind_network_combis_1, ind_network_combis_2] = find(triu(ones(nNetworks))==1); 
combis_networks = [ind_network_combis_1, ind_network_combis_2]; % all within and between network combinations

combis_states = nchoosek((1:nStates),2); %all state combinations


%% Compute relation between intelligence and cosine distance 
% distances between each within and between network combination for each
% combination of states
for n = 1:size(dist_cosine_join,1) % loop over network combinations

    
    % pairwise distance between all states for a network combination 
    dist_cosine_n = squeeze(dist_cosine_join(n,:,:));
    
    [rho_cosdist, p_cosdist] = partialcorr(g_score_sel,dist_cosine_n',confounds_sel,'type','spearman');
    rho_cosdist_all(:,n) = rho_cosdist;
    p_cosdist_all(:,n) = p_cosdist;    
    
end

% visualization
figure()
subplot(1,2,1)
h=imagesc(rho_cosdist_all);
set(h, 'AlphaData', ~isnan(rho_cosdist_all));
set(gca,'color',0.2*[1 1 1]);
caxis([-0.25 0.25]);
cMap = redblue;
colormap(flipud(cMap))  
bb = colorbar;
set(bb, 'ylim', [-0.25 0.15]);
%%% label single fields in matrix
% for n1=1:28
%     for n2 = 1:28
%         
% %         text(n2,sort_rho_corr_fcStatic_i(n1,n2), num2str(n1), 'HorizontalAlignment', 'center')
%     end
% end
xlabel('networks') 
ylabel('states') 
title('rho intelligence-reconfiguration')

% visualize which relations remain below a specific p-threshold

% maximal significant p-value after FDR 
[~, p_thresh_vis, ~, ~] = fdr_bh(p_cosdist_all(isfinite(p_cosdist_all)),alpha_FDR);

ww1 = size(rho_cosdist_all,1);
ww2 = size(rho_cosdist_all,2);

ind = find(p_cosdist_all(:) > p_thresh_vis);
rho_cosdist_all_tresh = rho_cosdist_all(:);
rho_cosdist_all_tresh(ind) = NaN;
rho_cosdist_all_tresh = reshape(rho_cosdist_all_tresh,ww1,ww2);

subplot(1,2,2)
[nr,nc] = size(rho_cosdist_all_tresh);
pcolor([rho_cosdist_all_tresh nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
bb = colorbar;
caxis([-0.25 0.25]);
set(bb, 'ylim', [-0.25 0.15]);
xlabel('networks') 
ylabel('states') 
title('rho intelligence-reconfiguration')

%% compute variance of correlations over task and network correlations 
% and test if there is a significant difference 
std_tasks = nanstd(rho_cosdist_all,0,2); % SD of rho between state combinations
std_networks = nanstd(rho_cosdist_all,0,1); % SD of rho between network combinations
figure()
hold on
histogram(std_tasks)
histogram(std_networks)

[pMWU,hMWU, statsU] = ranksum(std_networks,std_tasks); % test if SDs are significantly different
annotation('textbox',[.4 .5 .3 .3],'String',append('Wilcoxon rank sum test: W = ',num2str(statsU.ranksum),', p = ',num2str(round(pMWU,5))),'FitBoxToText','on');
xlabel('frequency') 
ylabel('variance') 
legend('SD over tasks','SD over networks')
title('histograms of variances in rho')
%% Relation between intelligence and total reconfiguration of within and between network reconfiguration (average over states) 
 
[rho_within_between_networks, p_within_between_networks] = partialcorr(g_score_sel,squeeze(nanmean(dist_cosine_join,2))',confounds_sel,'type','spearman');

% visualization
figure()
M_within_between_networks = zeros(nNetworks);
for e = 1:size(dist_cosine_join,1)
    M_within_between_networks(ind_network_combis_1(e),ind_network_combis_2(e)) = rho_within_between_networks(e); 
end

h = imagesc(M_within_between_networks);
line(repmat([1.5:1:nNetworks-0.5],2,1),repmat([0;nNetworks+1], 1, nNetworks-1), 'Color', 'black'); % vertical
line(repmat([0;nNetworks+1], 1, nNetworks-1), repmat([1.5:1:nNetworks-0.5],2,1), 'Color', 'black'); % horizontal
set(h, 'AlphaData', ~isnan(M_within_between_networks));
set(gca,'color',0.2*[1 1 1]);
caxis([-0.25 0.25]);
cMap = redblue;
colormap(flipud(cMap))  
bb = colorbar;
set(bb, 'ylim', [-0.25 0.12]);

% mark significance single fields in matrix
% maximal significant p-value after FDR 
[~, crit_p, ~, ~] = fdr_bh(p_within_between_networks(isfinite(p_within_between_networks)), alpha_FDR);

ind_sig = find(p_within_between_networks <= crit_p);
for e = 1: length(ind_sig)
      text(ind_network_combis_2(ind_sig(e)),ind_network_combis_1(ind_sig(e)), '*', 'HorizontalAlignment', 'center','fontsize',18)
end

xlabel('networks')
ylabel('networks') 
title('rho reconfiguration within/between networks - intelligence')

%% Total reconfiguration of each network across all states
for n=1:nNetworks
    
    ind = combis_networks(:,1) == n | combis_networks(:,2) == n;
    
    dist_cosine_join_reshaped = reshape(dist_cosine_join(ind,:,:),sum(ind)*size(dist_cosine_join,2),[]);
    [rho_total_networks(n), p_total_networks(n)] = partialcorr(g_score_sel,nanmean(dist_cosine_join_reshaped,1)',confounds_sel,'type','spearman');

end

% visualization
figure()

% maximal significant p-value after FDR 
[~, crit_p, ~, ~] = fdr_bh(p_total_networks(isfinite(p_total_networks)), alpha_FDR);  

labels = strings(length(p_total_networks));
for l = 1:length(p_total_networks)
    if p_total_networks(l) <= crit_p
        labels(l) = '*';
    else 
        labels(l) = '';
    end
end
plot_bar_annot(rho_total_networks, labels)
ylim([-0.25,0.05])
xlabel('networks') 
ylabel('rho')
title('rho network total reconfiguration - intelligence')
%% State specific relation between intelligence and network-wise reconfiguration

for c = 1:nStates
    
    ind_c = combis_states(:,1) == c | combis_states(:,2) == c; % indexes of state combinations of state c  
    
    if c>1
       ind_c(1:nStates-1)=0; % remove rest-task combination for task comparisons
    end
    
    for n = 1:nNetworks
        
        ind_n = combis_networks(:,1) == n | combis_networks(:,2) == n; % indexes of network combinations of network n
        
        dist_nc = dist_cosine_join(ind_n,ind_c,:); % all distances of specific state and network
        dist_nc_reshaped = reshape(dist_nc,sum(ind_n)*sum(ind_c),[]); % reshape vector
        mean_dist_nc = squeeze(nanmean(dist_nc_reshaped)); % average of distances
        [rho_dist_nc, p_dist_nc] = partialcorr(g_score_sel,mean_dist_nc',confounds_sel,'type','spearman');
        rho_dist_nc_all(n,c) = rho_dist_nc; % storing relations for each state and network
           
    end
end

% visualization
figure()
imagesc(rho_dist_nc_all')
caxis([-0.25 0.25]);
cMap = redblue;
colormap(flipud(cMap))  
bb = colorbar;
set(bb, 'ylim', [-0.25 0.12]);
xlabel('networks') 
ylabel('states') 
title('rho intelligence - network and state specific reconfiguration')