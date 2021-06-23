%%% Scope:  Script for selecting subjects according to completeness 
%%%         of fMRI data and cognitive scores, exclusion of high 
%%%         motion subjects
%%% Author: Jonas Thiele
%%% Date:   27.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% load parameter for analysis
load init_parameters


%% fMRI data avaiable?

participants_data = tdfread('participants.tsv'); % read subject IDs and data
subjects_IDs_all = participants_data.participant_id;
subjects_number_files_all = zeros(length(subjects_IDs_all),1); % vector for counting number of available fMRI files

for s=1:length(subjects_IDs_all) % loop over subjects to count number of fMRI files available
    
    subject_ID = subjects_IDs_all(s,:);
    folder_files_fMRI = fullfile('fmrip2mat_36pNS', subject_ID);
    filecount = 0;
    
    for sc = 1: length(scans)
        filepath = fullfile('fmrip2mat_36pNS', subject_ID, append(cell2mat(scans(sc)),'_out'),'output_makemat',filename_fMRI); 
        
        if isfile(filepath)
        % file exists
        filecount = filecount + 1;
        end
    end
    subjects_number_files_all(s) = filecount;
    
end
ind_fMRI_complete = subjects_number_files_all == length(scans); % indexes of subjetcs with complete fMRI

%% 
% Raven-score available?
raven_score = participants_data.raven_score;
ind_raven = nan(length(subjects_IDs_all),1); % indexes of subjects with raven score
for s=1:length(subjects_IDs_all)
    if strcmp(raven_score(s,:), append('n/a',blanks(1)))
        ind_raven(s) = 0;
    else
        ind_raven(s) = 1;
    end
end


% confounds available?

age = participants_data.age;
ind_age = nan(length(subjects_IDs_all),1); % indexes of subjects with age value
for s=1:length(subjects_IDs_all)
    if strcmp(age(s,:), append('n/a',blanks(2)))
        ind_age(s) = 0;
    else
        ind_age(s) = 1;
    end
end

sex = participants_data.sex;
ind_sex = nan(length(subjects_IDs_all),1); % indexes of subjects with sex value
for s=1:length(subjects_IDs_all)
    if strcmp(sex(s,:), 'n/a')
        ind_sex(s) = 0;
    else
        ind_sex(s) = 1;
    end
end

handedness = participants_data.handedness;
ind_handedness = nan(length(subjects_IDs_all),1); % indexes for subjects with handedness value
for s=1:length(subjects_IDs_all)
    if strcmp(handedness(s,:), append('n/a',blanks(9)))
        ind_handedness(s) = 0;
    else
        ind_handedness(s) = 1;
    end
end
    
  
% motion exclusion

subject_motion_all = nan(length(subjects_IDs_all),length(scans)*3);
subject_high_motion_all = nan(length(subjects_IDs_all),length(scans)*3);
threshSmallSpike = 0.25;
threshMeanFD = 0.2;
threshLargeSpike = 5;
threshPercent = 0.2;

for s=1:length(subjects_IDs_all)
    
    subject_ID = subjects_IDs_all(s,:);
    folder_files_regressors = fullfile('regressors',subject_ID,'func'); % folder with head motion (framewise displacment) data
    
    for sc = 1: length(scans)
        filename = append(subject_ID,'_',cell2mat(scans(sc)),'_desc-confounds_regressors.tsv'); 
        filepath = fullfile(folder_files_regressors,filename); 
        
        if isfile(filepath)
            FD = tdfread(filepath).framewise_displacement; % get head motion for subjects s of scan sc
            FD = str2num(FD(2:end,:));
           
            numberTimeSteps = length(FD); % length of scan
            meanFD = mean(FD); % mean framewise displacement
            numSmallSpikes = sum(FD > threshSmallSpike); % number of small spikes
            numLargeSpikes = sum(FD > threshLargeSpike); % number of large spikes
            PercSmallSpikes = numSmallSpikes/numberTimeSteps; % percentage of small spikes 
            
            % store mean FD, small spikes, large spikes in
            % subject_motion_all
            subject_motion_all(s,(sc-1)*3+1) = meanFD;
            subject_motion_all(s,(sc-1)*3+2) = PercSmallSpikes;
            subject_motion_all(s,(sc-1)*3+3) = numLargeSpikes;
            
            % check if FD parameters exceed thresholds and store information in subject_high_motion_all 
            if meanFD >= threshMeanFD
                subject_high_motion_all(s,(sc-1)*3+1) = 1;
            else
                subject_high_motion_all(s,(sc-1)*3+1) = 0;
            end
            
            if PercSmallSpikes >= threshPercent
                subject_high_motion_all(s,(sc-1)*3+2) = 1;
            else
                subject_high_motion_all(s,(sc-1)*3+2) = 0;
            end
            if numLargeSpikes > 0
                
                subject_high_motion_all(s,(sc-1)*3+3) = 1;
            else
                subject_high_motion_all(s,(sc-1)*3+3) = 0;
            end
                      
        end
    end

    
end    


ind_motion = sum(subject_high_motion_all,2) == 0; % indexes for fullfilment of motion criteria

% index of subjects for which all begavioral and fMRI data is avaiable 
% and motion is within criteria
ind_all_ok = ind_fMRI_complete.*ind_age.*ind_sex.*ind_handedness.*ind_raven.*ind_motion; 

%% get IDs, raven-score, and confounds for subjects that meet all criteria 

subject_IDs_sel = [];
cnt = 0;
for s=1:length(subjects_IDs_all)
   
    if ind_all_ok(s) == 1 % subjects that meet all criteria
       cnt = cnt + 1; 
       subject_IDs_sel = [subject_IDs_sel; {subjects_IDs_all(s,:)}]; 
       raven_score_sel(cnt) = str2num(raven_score(s,:));
       
       age_sel(cnt) = str2num(age(s,:)); % age value
       
        % M = 1, F = 0;
        if strcmp(strtrim(sex(s,:)),'M') 
            sex_s = 1;
        elseif strcmp(strtrim(sex(s,:)),'F') 
            sex_s = 0;
        end
        sex_sel(cnt) = sex_s; % gender
       
        % -1 = left, 1=right, 0=ambidextrous
        if strcmp(strtrim(handedness(s,:)),'right') 
            hand_s = 1;
        elseif strcmp(strtrim(handedness(s,:)),'left') 
            hand_s = -1;
        elseif strcmp(strtrim(handedness(s,:)),'ambidextrous')
            hand_s = 0; 
        end
        handedness_sel(cnt) = hand_s; % handedness
        
        % individual's average of mean framewise displacement across all scans
        meanFD_allScans_sel(cnt) = mean(squeeze(subject_motion_all(s,1:3:end))); 
        
        % individual's average of percentage of small spikes across all scans
        mean_perc_small_spikes(cnt) = mean(squeeze(subject_motion_all(s,2:3:end)));
        
    end
    
end

%% make a table with all variables of selected subjects
data_beh_sel = table(subject_IDs_sel,raven_score_sel',age_sel',sex_sel',...
    handedness_sel',meanFD_allScans_sel',mean_perc_small_spikes');

% translate the variables to a suitable format 
names_variables = {'Subject','raven_score','age','sex','handedness','meanFD','perc_small_spikes'};
data_beh_sel.Properties.VariableNames = names_variables;

writetable(data_beh_sel, 'data_beh_sel_subjects.csv');

% save confounds and intelligence score also separately
confounds_sel = [age_sel',sex_sel',...
    handedness_sel',meanFD_allScans_sel',mean_perc_small_spikes'];
g_score_sel = raven_score_sel';

save('confounds_sel.mat', 'confounds_sel')
save('g_score_sel.mat', 'g_score_sel')