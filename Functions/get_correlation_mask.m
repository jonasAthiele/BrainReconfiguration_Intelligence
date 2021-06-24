%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coded by: Jonas Thiele
%%% Date: 24.06.2021
%%% Platform: MATLAB
%%% Purpose: get mask of edges correlating positively or negatively (as per sign)  
%%%          with the g-score (below specific p-threshold) 
%%%  	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function corr_mask = get_correlation_mask(FC, g_score, confounds, p_tresh, nScan, sign)
	
% Input Arguments
%__________________________________________________________________________
%
%   FC          -- correlation matrices of the sample 
%   g_score     -- intelligence scores of the sample
%   confounds   -- confounds of the sample
%   p_tresh     -- treshhold for filtering
%   nScan       -- number of scan
%   sign        -- pos (1) or neg (-1) for direction of correlation 
%
% Output Arguments
%__________________________________________________________________________
%
%  mask_corr    -- edges correlating positively or negatively (as per sign) 
%                  with the g-score (below p_tresh) 
%                  

%%
    [nNodes, ~, nSubjects, ~] = size(FC); % number of nodes and subjects
    
    FC_state = squeeze(FC(:,:,:,nScan)); % FC of specific state
    FC_state = reshape(FC_state, nNodes*nNodes, nSubjects); % reshape FC
    
    % correlating connection strengths with intelligence
    [rho_fc, p_fc] = partialcorr(g_score,FC_state',confounds,'type','spearman');
    
    % find indexes of connections correlating with intelligence with p
    % smaller threshhold
    index_p1 = abs(p_fc) < p_tresh;
    
    % only keep a specific direction of correlations (positive or negative)
    if sign == 1
        index_p2 = rho_fc > 0;
    elseif sign == -1
        index_p2 = rho_fc < 0;
    else 
        error('sign must be either 1 or -1')
    end
    
    % mask with connections correlating positively or negatively 
    % (dependend on sign) and with p smaller a p_thresh with intelligence 
    corr_mask = index_p1.*index_p2;
    corr_mask = corr_mask';
    
end