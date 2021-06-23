

%%% Scope:  Create samples  
%%%         with absence of family relations between subjects
%%%         and an equal distribution of intelligence scores via 
%%%         stratified folds 
%%%         + save corresponding confounds and intelligence scores
%%% Author: Jonas Thiele
%%% Date:   16.06.2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%% 1: Find families
% load parameters of analysis
load init_parameters
% read behavioral data of subject selection
beh_sel = readtable('data_beh_sel_subjects.csv');
% read restricted data of all subjects
res_all = readtable('tables_HCP/RESTRICTED_jonasthiele_11_1_2020_8_46_1.csv');

% choose restricted data from relevant subjects only
subs_sel = beh_sel.Subject;
subs_all = res_all.Subject;
ind_sel = ismember(subs_all, subs_sel);
res_sel = res_all(ind_sel,:);

% mother and father IDs of all selected subjects
ID_vec_mom = res_sel.Mother_ID;
ID_vec_dad = res_sel.Father_ID;

% look for subjects with similar mother or father ID and join them to
% families
cnt = 1;
family_vec = [];
for s=1:length(res_sel.Subject)
    

    ID_mom_s = res_sel.Mother_ID(s);
    ID_dad_s = res_sel.Father_ID(s);
    
    ind_moms = find(res_sel.Mother_ID == ID_mom_s);
    ind_dads = find(res_sel.Father_ID == ID_dad_s);
    
    family_temp = unique([ind_moms; ind_dads]);
    
    % check if family already exists
    repeat = 0;
    for c = 1:cnt-1
    
        fam_old_members = family_vec(find(family_vec(:,2) == c),1);
     
            
        if isequal(sort(fam_old_members),sort(family_temp))
            repeat = 1;
        end
        
    end
        
    if repeat == 0
        family_vec = [family_vec; [family_temp, ones(length(family_temp),1)*cnt]];
        cnt = cnt + 1;
    end
end

%% Looking for doubling subjects --> subjects can occur in multible families
subjects_recollect = [];
subjects_doubling = [];
for s=1:length(family_vec)
    
   if ismember(family_vec(s,1), subjects_recollect)
      subjects_doubling = [subjects_doubling;family_vec(s,1)];
   end
    
   subjects_recollect = [subjects_recollect; family_vec(s,1)];
    
end
subjects_doubling = unique(subjects_doubling);

%% Joining families 
% all subjects that are somehow related (at least one same subject multiple families) are joined into one family
% e.g. families [sub1,sub2] and [sub1,sub3] are joined to [sub1, sub2, sub3]

family_vec_joined = family_vec;
for s = 1: length(subjects_doubling)
   
    ind = find(family_vec_joined(:,1) == subjects_doubling(s));
    
     
    if length(ind) > 1
        run = 1;
        % get the whole chain --> all subjetcs somehow connected should be
        % found
        elements = subjects_doubling(s); % in elements all subjects that are somehow connected are stored
        fam_no = []; % the numbers of the families of all connected subjects are stored here
        elements_old = []; % for comparing if new subjects were added to elements
        while run == 1
            
            for e=1:length(elements)
               
                ind_e = find(family_vec(:,1) == elements(e));
                
                fam_no = unique([fam_no; family_vec(ind_e,2)]); % the families of subjects in elements 
                
            end
            
            for n=1:length(fam_no)
                
                elements = unique([elements; family_vec(find(family_vec(:,2)==fam_no(n)),1)]);
            end
             
            if isequal(elements,elements_old)
                run = 0;
            end
            
            elements_old = elements;
        end
    
        
        % families to join were collected in fam_no, subjects to join are in
        % elements
        fam_no_joined = fam_no(1); % name of joined family
        
        % delete all families that are joined
        for f = 1:length(fam_no)
           family_vec_joined(find(family_vec_joined(:,2)==fam_no(f)),:) = []; 
        end
        
        % create new family
        family_vec_joined = [family_vec_joined; [elements, ones(length(elements),1)*fam_no_joined]];
        
        % sorting families according to their number
        [vec_sorted, vec_order] = sort(family_vec_joined(:,2));
        family_vec_joined = family_vec_joined(vec_order,:);
    end
end
% renumbering families
fam_no_unique = unique(family_vec_joined(:,2));
for f = 1:length(fam_no_unique)
    family_vec_joined(family_vec_joined(:,2) == fam_no_unique(f),2) = f;
end
    
%% adding IDs to each subject to better compare 
for s = 1:length(family_vec_joined)
    
    mom_vec_joined(s) = ID_vec_mom(family_vec_joined(s,1));
    dad_vec_joined(s) = ID_vec_dad(family_vec_joined(s,1));
end
family_vec_joined = [family_vec_joined, mom_vec_joined', dad_vec_joined'];


%% Calculating g-values for each family
% get g-scores of selected subjects
g_score_1186 = readtable('gfac_bi_explore_1186Subjects.csv').x;
subjects_1186 = readtable('subject_IDs_g_score_1186Subjects.csv').Subject;
ind_sel = ismember(subjects_1186, res_sel.Subject);
g_score_sel = g_score_1186(ind_sel);

% add g-value to each subject
for s=1:length(family_vec_joined)
    g_vec_joined(s)=g_score_sel(family_vec_joined(s,1));
end

% averaging the g-values within each family to get an family average g
g_average_vec_joined = [];
for f=1:max(family_vec_joined(:,2))
    ind = family_vec_joined(:,2) == f;
    g_average_vec_joined = [g_average_vec_joined; ones(sum(ind),1)*mean(g_vec_joined(ind))];
end

family_vec_joined = [family_vec_joined, g_vec_joined', g_average_vec_joined];

%% Sort families according to their g-score
families_g = [];
for f=1:max(family_vec_joined(:,2))
    ind = find(family_vec_joined(:,2) == f);
    families_g = [families_g; [family_vec_joined(ind(1),2),family_vec_joined(ind(1),6)]]; 
    
end
[vec_sorted, vec_order] = sort(families_g(:,2));
families_g = families_g(vec_order,:);

% split into samples
% get all families per sample
for s=1:nSamples
    %families_samples(:,f) = families_g(f:nSamples:end,1);
    families_samples(s).members = families_g(s:nSamples:end,1);
end

% get all subjects per sample
for s=1:nSamples
    subjects_s = [];
    for f=1:length(families_samples(s).members)
        
        ind = find(family_vec_joined(:,2) == families_samples(s).members(f)); 
        
        subjects_f = family_vec_joined(ind,1); 
        subjects_s = [subjects_s; subjects_f];
        
        
        
    end
    
    % store subjects, IDs of subjects and g-scores of respective sample in a struct
    subjects_s = sort(subjects_s);
    subjects_samples(s).members = subjects_s;
    subjects_samples(s).IDs = subs_sel(subjects_s);
    subjects_samples(s).g_score = g_score_sel(subjects_s);
    
end


%% Get confounds for selected subjects

FD_mean_data = readtable('dataMeanFD_allfMRIsets_948Subjects.csv');
spikes_data = readtable('dataSmallSpikes_allfMRIsets_948Subjects.csv');

ind_sel = ismember(FD_mean_data.Subjects, res_sel.Subject);
FD_mean_data_sel = FD_mean_data(ind_sel, :);

ind_sel = ismember(spikes_data.Subjects, res_sel.Subject);
spikes_data_sel = spikes_data(ind_sel, :);

age_sel = res_sel.Age_in_Yrs;
hand_sel = res_sel.Handedness;
gender_sel = cell2mat(beh_sel.Gender) == 'M';

mean_meanFD = table2array(FD_mean_data_sel); 
mean_meanFD = mean(mean_meanFD(:,2:end),2); %average of the mean framewise displacements across all scans

mean_spikes = table2array(spikes_data_sel);
mean_spikes = mean(mean_spikes(:,2:end),2); %average of percentage of small spikes across all scans


confounds_sel = horzcat(age_sel, hand_sel, gender_sel, mean_meanFD, mean_spikes); 

% add confounds to the struct
for s = 1:length(subjects_samples)

    subjects_samples(s).confounds = confounds_sel(subjects_samples(s).members,:);
    
end


%% save data of samples and additional confounds and g-scores 
save('subjects_samples.mat','subjects_samples') 
save('confounds_sel.mat', 'confounds_sel')
save('g_score_sel.mat','g_score_sel')
%% Distribution of g-scores in the samples

% figure()
% for k=1:nSamples
%    
%     subplot(2,5,k)
%     
%     histogram(subjects_samples(k).g_score,10)
%     
% end