

%% Function for plottting matrices with correlations between reconfiguration and intelligence
% written by Jonas Thiele
% date: 24.06.2021

function plot_rho_p_fdrCorrected(rho, p, alpha_FDR, nCondis, range, varName)

%% Input:
%  rho...array of correlations
%  p...array of p-values of correlations rho
%  alpha_FDR...significance level (e.g. 0.05)
%  nCondis...number of tasks that are compared
%  varName...Name of reconfiguration measure (e.g. 'cosine distance')
%  range...range of values shown in the colorbar of the plot
%% Output: 
%  pFDR: p...largest significant p-value (result of FDR)

%%
combis_tasks = nchoosek(1:nCondis,2);


rho_M=NaN(nCondis);

for c = 1:length(combis_tasks)
    rho_M(combis_tasks(c,1),combis_tasks(c,2))= rho(c);
end

h=imagesc(rho_M);
line(repmat([1.5:1:nCondis-0.5],2,1),repmat([0;nCondis+1], 1, nCondis-1), 'Color', 'black'); % vertical
line(repmat([0;nCondis+1], 1, nCondis-1), repmat([1.5:1:nCondis-0.5],2,1), 'Color', 'black'); % horizontal
ax = gca;
ax.TickLength = [0 0];
caxis([-max(abs(range)),max(abs(range))]); % centering white color
set(h, 'AlphaData', ~isnan(rho_M));
cMap = redblue;
colormap(flipud(cMap))  
bb = colorbar;
set(bb, 'ylim', range)
[~, crit_p, ~, ~] = fdr_bh(p(isfinite(p)), alpha_FDR); % smallest significant p-value after FDR correction 
for c = 1: 1:length(combis_tasks)
    if p(c) <= crit_p 
       text(combis_tasks(c,2),combis_tasks(c,1), '*', 'HorizontalAlignment', 'center','fontsize',18)
    end
end
str_title = append('rho: ', varName);
xlabel('states') 
ylabel('states') 
title(str_title)


