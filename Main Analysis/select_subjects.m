
%%% Scope:  Script for the selection of subjects according to data 
%%%         availability and motion exclusion criteria
%%% Author: Jonas Thiele
%%% Date:   16.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% load parameters of analysis
load init_parameters

% behavioral data provided by the HCP
dataBeh = readtable('tables_HCP/unrestricted_jonasthiele_11_1_2020_8_39_10.csv');


%% Level 1 exclusion: exclude subjects without PMAT + 11 other task scores

dataSel=dataBeh((~isnan(dataBeh.PicVocab_Unadj)) &...
                (~isnan(dataBeh.ReadEng_Unadj)) &...
                (~isnan(dataBeh.PicSeq_Unadj)) &...
                (~isnan(dataBeh.Flanker_Unadj)) &...
                (~isnan(dataBeh.CardSort_Unadj)) &...
                (~isnan(dataBeh.ProcSpeed_Unadj)) &...
                (~isnan(dataBeh.PMAT24_A_CR)) &...
                (~isnan(dataBeh.VSPLOT_TC)) &...
                (~isnan(dataBeh.IWRD_TOT)) &...
                (~isnan(dataBeh.ListSort_Unadj)) &...
                (~isnan(dataBeh.DDisc_AUC_200)) &...
                (~isnan(dataBeh.DDisc_AUC_40K)) &...
                (~isnan(dataBeh.SCPT_TP))&...
                (~isnan(dataBeh.SCPT_TN))&...  
                (~isnan(dataBeh.SCPT_FP))&...  
                (~isnan(dataBeh.SCPT_FN))&...
                (~isnan(dataBeh.SCPT_TPRT)),:);


dataSel=dataSel(strcmp(dataSel.DelDisc_Compl,'true') &...
                 strcmp(dataSel.SCPT_Compl,'true') &...
                 strcmp(dataSel.IWRD_Compl,'true') &...
                 strcmp(dataSel.VSPLOT_Compl,'true'),:);               
            
%% Level 2 exclusion: exclude subjects with MMSE score <= 26
dataSel=dataSel(strcmp(dataSel.MMSE_Compl,'true') & (dataSel.MMSE_Score >=26),:);  

%% Compute combined SCPT and combined DDisc measures

SCPT_TP = normalize(dataSel.SCPT_TP,'range',[0,1]);
SCPT_TN = normalize(dataSel.SCPT_TN,'range',[0,1]);
SCPT_FP = normalize(dataSel.SCPT_FP,'range',[0,1]);
SCPT_FN = normalize(dataSel.SCPT_FN,'range',[0,1]);
SCPT_TPRT = normalize(dataSel.SCPT_TPRT,'range',[1,2]);
DDisc_AUC_200 = normalize(dataSel.DDisc_AUC_200,'range',[0,1]);
DDisc_AUC_40K = normalize(dataSel.DDisc_AUC_40K,'range',[0,1]);

SCPT_Eff = (SCPT_TP + SCPT_TN)./ ((SCPT_TP + SCPT_TN + SCPT_FN + SCPT_FP).*SCPT_TPRT);
DDisc = DDisc_AUC_200 + DDisc_AUC_40K;


%% Put all 12 cognitive scores with combined SCPT and DDisc scores in a matrix

dataBeh12 = [dataSel.PicVocab_Unadj,dataSel.ReadEng_Unadj,dataSel.PicSeq_Unadj,...
    dataSel.Flanker_Unadj,dataSel.CardSort_Unadj,dataSel.ProcSpeed_Unadj,dataSel.PMAT24_A_CR,...
    dataSel.VSPLOT_TC,dataSel.IWRD_TOT,dataSel.ListSort_Unadj,SCPT_Eff,DDisc];

dataBeh12 = normalize(dataBeh12);

dataBeh12_names = {'PicVocab_Unadj','ReadEng_Unadj','PicSeq_Unadj','Flanker_Unadj',...
                             'CardSort_Unadj','ProcSpeed_Unadj','PMAT24_A_CR','VSPLOT_TC','IWRD_TOT',...
                             'ListSort_Unadj','SCPT_Eff','DDisc'};
                         
T_dataBeh12 = array2table(dataBeh12);
T_dataBeh12.Properties.VariableNames = dataBeh12_names;

T_subject_ID_g_score = array2table(dataSel.Subject);
T_subject_ID_g_score.Properties.VariableNames = {'Subject'};
writetable(T_subject_ID_g_score, 'subject_IDs_g_score_1186Subjects.csv')

% g-scores are computed from these cognitive scores in R using bi-factor analysis model 

%% Level 3 exclusion: take only data with all rest and task fMRI

dataSel=dataSel((strcmp(dataSel.x3T_Full_Task_fMRI, 'true')) & (dataSel.x3T_RS_fMRI_PctCompl == 100.0),:);
dataSel=dataSel(dataSel.Subject~=668361,:);  % exclude subject manually because tfRMI_WM_RL missing
%% Level 4 exclusion: subjects with high motion

% thresholds for highMotion 
threshSmallSpike = 0.25; % framewise displacement over thresh is a small spike
threshMeanFD = 0.2; % maximal allowed mean freamwise displacement
threshLargeSpike = 5; % framewise displacement over thresh is a large spike
threshPercent = 0.2; % maximal allowed percentage of small spikes in a scan

% path of Movement_RelativeRMS_files folder(Files with movements during scan of the subjects)
folder = 'Movement_RelativeRMS_files';

nSubjects = length(dataSel.Subject); % number of subjects
nScans = length(scans); % number of scans

% initialize vectors for storing motion data
dataRMSallSubj = [];
highMotionAllSubj = [];

for s = 1:nSubjects
    
    ID_subject = dataSel.Subject(s); % ID of subjects
    dataRMSSubj=[];
    highMotionSubj=[];
    for sc = 1:nScans
        
        % path of the RMS file for the specific subject for the specific task
        pathRMSdata=fullfile(folder, string(ID_subject),'MNINonLinear','Results',scans(sc));
        
        % open, read the .txt file and save it in an array 
        file = fullfile(pathRMSdata,'Movement_RelativeRMS.txt');
        fileID = fopen(file,'r');
        formatSpec = '%f';
        relRMS =  fscanf(fileID,formatSpec);
        fclose(fileID);
        
        % compute relevant RMS measures
        numberTimeSteps = length(relRMS);
        meanFD = mean(relRMS);
        numSmallSpikes = sum(relRMS > threshSmallSpike);
        numLargeSpikes = sum(relRMS > threshLargeSpike);
        PercSmallSpikes = numSmallSpikes/numberTimeSteps;
        
        % store RMS values extracted from RMS file in dataRMSSubj array
        dataRMSSubj = [dataRMSSubj, [meanFD,numSmallSpikes,numLargeSpikes]];
        
        % check if values above treshhold and save results in highMotionSub array
        highMotionSubj = [highMotionSubj, sum(meanFD>=threshMeanFD)];
        highMotionSubj = [highMotionSubj, sum(PercSmallSpikes>=threshPercent)];
        highMotionSubj = [highMotionSubj, sum(numLargeSpikes>0)];
        
    end
    
    % dataRMSSubj = [ID_subject, dataRMSSubj]; % add ID of subject into array
    dataRMSallSubj = [dataRMSallSubj;dataRMSSubj];
    highMotionAllSubj = [highMotionAllSubj; highMotionSubj];  
   
end

% make tables with mean framwise displacement and spikes data
T_dataMeanFD = array2table([dataSel.Subject, dataRMSallSubj(:, 1:3:end)]);
T_dataMeanFD.Properties.VariableNames = ['Subjects',scans];

T_dataSmallSpikes = array2table([dataSel.Subject, dataRMSallSubj(:, 2:3:end)]);
T_dataSmallSpikes.Properties.VariableNames = ['Subjects',scans];

% exclude all subjects with high motion = final sample
data_beh_sel = dataSel(sum(highMotionAllSubj,2)==0,:);

% save behavioral and motion data 
writetable(data_beh_sel, 'data_beh_sel_subjects.csv');
writetable(T_dataBeh12, 'cogScores_1186Subjetcs.csv');

writetable(T_dataMeanFD, 'dataMeanFD_allfMRIsets_948Subjects.csv');
writetable(T_dataSmallSpikes, 'dataSmallSpikes_allfMRIsets_948Subjects.csv');
