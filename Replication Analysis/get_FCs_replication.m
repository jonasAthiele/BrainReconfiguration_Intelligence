
%%% Scope:  Reading preprocessed BOLD-signals, calculating FCs 
%           (functional connectivities), joining FCs
%%% Author: Jonas Thiele
%%% Date:   17.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
% load parameters for analysis
load init_parameters
 
% load subjects
behavioralData = readtable('data_beh_sel_subjects.csv');
subjects = behavioralData.Subject;

nSubjects = length(subjects);

for sc = 1:nStates
    
    for s=1:nSubjects % numberSubjects 

        idSubject = cell2mat(subjects(s));
        
        % read fMRI data
        filepath = fullfile('fmrip2mat_36pNS', idSubject, append(cell2mat(scans(sc)),'_out'),'output_makemat',filename_fMRI); 
        
        file = h5read(filepath, '/timeseries');
      
        bold = squeeze(file); % region x neural activity

        % remove background nodes 
        bold(background_regions,:)=[];
        bold=bold(1:nNodes,:);

        FCi = corr(bold'); % correlating node signals

        FCi = (FCi+FCi')./2; % symmetrize matrix
        FCi(1:size(FCi,1)+1:end) = 0; % set diagonal elements to zero
        FCi = fisherZTransform(FCi); % fisher Z transform all correlations

        FCStatic_combined(:,:,s,sc) = FCi; % named as combined for applicability on codes from HCP sample

    end
    
end

save('FCStatic_combined.mat', 'FCStatic_combined')



