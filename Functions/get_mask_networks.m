%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coded by: Jonas Thiele 
%%% Date: 24.06.2020
%%% Platform: MATLAB
%%% Purpose: Compute masks of inter network and pairwise intra network connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mask_networks = get_mask_networks(node_assignments) 	% Input is the correlation matrix C

% Input Arguments
%__________________________________________________________________________
%
% nodes_assignments -- 1D-vector of network labels for each node  
%
% Output Arguments
%__________________________________________________________________________
%
%   mask_networks -- masks of nodes for each within and between network
%                    combination (nodes x nodes x network combinations)



nNetworks = max(node_assignments); % number of networks
nNodes = size(node_assignments,1); % number of nodes

% indexes of all within and between network combinations
[x_ind, y_ind] = find(triu(ones(nNetworks))==1); 

% mask for storing within and between network connections
mask_networks = zeros(nNodes,nNodes, length(x_ind));

for n = 1:length(x_ind)
        
        
        nodes1 = find(node_assignments == x_ind(n)); % nodes network 1
        nodes2 = find(node_assignments == y_ind(n)); % nodes network 2
        % for within networks nodes1 = nodes 2

        nodes = unique([nodes1;nodes2]); % sort nodes of both networks
        combis_nodes = nchoosek(nodes,2); % all combinations of nodes   
        
        % set all node combionations of network combination to "1"
        for c=1:size(combis_nodes,1)
            mask_networks(combis_nodes(c,1),combis_nodes(c,2),n) = 1; 

        end
        % delete within network nodes for between network comparisons
        if x_ind(n) ~= y_ind(n)
            
            nodes = sort(nodes1);
            combis_nodes = nchoosek(nodes,2);
            
            for c=1:size(combis_nodes,1)
                mask_networks(combis_nodes(c,1),combis_nodes(c,2),n) = 0; 
            end
            nodes = sort(nodes2);
            combis_nodes = nchoosek(nodes,2);
            
            for c=1:size(combis_nodes,1)
                mask_networks(combis_nodes(c,1),combis_nodes(c,2),n) = 0; 
            end
            
        end
        mask_networks(:,:,n) = mask_networks(:,:,n) - diag(diag(mask_networks(:,:,n)));
end

