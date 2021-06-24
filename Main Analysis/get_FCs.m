
%%% Scope:  Reading preprocessed BOLD-signals, calculating FCs
%%%         (functional connectivities), joining FCs
%%% Author: Jonas Thiele
%%% Date:   24.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% load parameter for analysis
load init_parameters

% read behavioral data of subjects
behavioralData = readtable('data_beh_sel_subjects.csv');
subjects = behavioralData.Subject; % subject IDs

nSubjects = length(subjects); % number of subjects


for sc = 1:length(scans) % loop over all scans
    
    % choose folder in which files are stored - here deiffernt for rest and
    % tasks data
    if ismember(sc, inds_rest)
       folderData = folderRest;
    else 
       folderData = folderTask; 
    end
    
    for s=1:nSubjects % numberSubjects 

        idSubject = subjects(s);
        %%
        %--1 Read fMRI data

        filepath = fullfile(folderData, num2str(idSubject), string(scans(sc)), filename_nii);

        file = niftiread(filepath);
        bold = squeeze(file); % neural activity time x region

        bold=bold(:,1:nNodes); % select cortical nodes only

        FCi = corr(bold); % corellating time series of nodes

        FCi = (FCi+FCi')./2; % symmetrize matrix
        FCi(1:size(FCi,1)+1:end) = 0; % set diagonal elements to zero
        FCi = fisherZTransform(FCi); % fisher z-transform all correlations

        FCStatic(:,:,sc,s) = FCi;

    end
    
end

%% Averaging FCs over cognitive states (resting state, different tasks)
for sc = 1:nStates
      
    FCStatic_combined(:,:,:,sc) = squeeze(mean(FCStatic(:,:,mask_states==sc,:),3));
                 
end
%%  
clear FCStatic

save('FCStatic_combined.mat', 'FCStatic_combined')
