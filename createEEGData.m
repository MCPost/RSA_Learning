%% Import and Create Data Set

function Data_EEG = createEEGData(cfg)

%% cfg Struct

Name = cfg.Name;
chan_idx = cfg.chan_idx; 
Art_corr = cfg.Art_corr; 
BL_corr = cfg.BL_corr; 
BL_wind = cfg.BL_wind;


%% Load Data
fprintf('\n')
fprintf('Load Data for Subject %s created!', Name)
fprintf('\n')
Preproc_Data_1 = load(['Prepro_2_encoding_full_',Name,'.mat'], ['r_1_',Name]);
Preproc_Data_1 = Preproc_Data_1.(cell2mat(fieldnames(Preproc_Data_1)));
Preproc_Data_2 = load(['Prepro_2_encoding_full_',Name,'.mat'], ['r_2_',Name]);
Preproc_Data_2 = Preproc_Data_2.(cell2mat(fieldnames(Preproc_Data_2)));

logfile = importdata(['Logs/resultfile_',Name,'.csv']);

logfile.textdata(2:end,find(strcmp(logfile.textdata(find(strcmp(logfile.textdata(:,2),'retrieval'),1,'first'),:),''),1,'first')+(0:size(logfile.data,2)-1)) = num2cell(logfile.data);
logfile = logfile.textdata;
for i = 1:size(logfile,2)
    if(~isempty(str2double(logfile{2,i})) && ~isnan(str2double(logfile{2,i})))
        logfile(2:end,i) = num2cell(str2double(logfile(2:end,i)));
    end
end
clear i


%% Get Encoding from logfile and Categories of Objects

names = {'block_number','trial_encoding','trial_total_enc','object_name','object_number','object_percep',... 
         'object_category','RT_encoding','cue_retrieval','OnsetEncodingBlockStart','WordOnsetEncoding',... 
         'ObjectContextOnset','ResponseEncodingOnset'};

Encoding_Data = logfile(strcmp(logfile(:,2),'encoding'), cell2mat(cellfun(@(x) ismember(x,names), strtrim(logfile(1,:)),'UniformOutput', 0)));
% block_number, trial_encoding, trial_total_enc, object_name,
% object_number, object_percep, object_category, RT_encoding,
% cue_retrieval, OnsetEncodingBlockStart, WordOnsetEncoding, 
% ObjectContextOnset, ResponseEncodingOnset, PerceptualDim, SemanticDim,
% Category

Encoding_Data(cell2mat(Encoding_Data(:,6)) == 1,14) = {'Drawing'};
Encoding_Data(cell2mat(Encoding_Data(:,6)) == 2,14) = {'Picture'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) < 5,15) = {'Animate'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) > 4,15) = {'Inanimate'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) < 5,16) = {1};
Encoding_Data(cell2mat(Encoding_Data(:,7)) > 4,16) = {2};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 1,17) = {'insect'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 2,17) = {'bird'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 3,17) = {'mammal'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 4,17) = {'sea'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 5,17) = {'fruit'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 6,17) = {'electronic'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 7,17) = {'veggie'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 8,17) = {'clothe'};
cfg.Encoding_Data = Encoding_Data;

%% Marry Behavioral and EEG Data
fprintf('\n')
fprintf('Marry Behavioral and EEG Data for Subject %s!', Name)
fprintf('\n')
TrialInfo = cell(size(Encoding_Data,1),12);
Preproc_Data_2.trialinfo = num2cell(Preproc_Data_2.trialinfo);
draw_idx = find(cell2mat(Encoding_Data(:,6)) == 1); cnt = 1;
for i = 1:length(draw_idx)
    if(ismember(draw_idx(i), cell2mat(Preproc_Data_2.trialinfo(:,2))))
        TrialInfo(i,:) = [Preproc_Data_2.trialinfo(cnt,1:4) Encoding_Data(draw_idx(i),[6 14 16 15 7 17 5 4])];
        cnt = cnt + 1;
    else
        TrialInfo(i,:) = [{0}, {draw_idx(i)} Encoding_Data(draw_idx(i),[6 7 6 14 16 15 7 17 5 4])];
    end
end
Preproc_Data_1.trialinfo = num2cell(Preproc_Data_1.trialinfo);
pic_idx = find(cell2mat(Encoding_Data(:,6)) == 2); cnt = 1;
for i = 1:length(pic_idx)
    if(ismember(pic_idx(i), cell2mat(Preproc_Data_1.trialinfo(:,2))))
        TrialInfo(64+i,:) = [Preproc_Data_1.trialinfo(cnt,1:4) Encoding_Data(pic_idx(i),[6 14 16 15 7 17 5 4])];
        cnt = cnt + 1;
    else
        TrialInfo(64+i,:) = [{0}, {pic_idx(i)} Encoding_Data(pic_idx(i),[6 7 6 14 16 15 7 17 5 4])];
    end
end

% Get Index of NaN trials
Nan_trial_idx = find(cell2mat(TrialInfo(:,1)) == 0);


%% Create EEG Data Matrix
[~,sort_idx] = sortrows(cell2mat(TrialInfo(:,[5 7 9])));

cfg.TrialInfo = TrialInfo(sort_idx,:);
cfg.TimeVec1024 = Preproc_Data_1.time;

Data = [Preproc_Data_2.trial; Preproc_Data_1.trial];

if(Art_corr)
    fprintf('\n')
    fprintf('Artefact Rejection for Subject %s!', Name)
    fprintf('\n')
    load ArtifactRej
    for i = 1:size(trial2rm.(Name),1)
        if(~strcmp(trial2rm.(Name){i,2},'all'))
            Data(trial2rm.(Name){i,2},trial2rm.(Name){i,1},:) = NaN;
        end
    end
    cfg.chan_idx_usable = chan_idx(~ismember(chan_idx, cell2mat(trial2rm.(Name)(strcmp(trial2rm.(Name)(:,2),'all'),1))));
end

% Baseline Correction
[~,BL_wind_idx] = min(abs(bsxfun(@minus, [cfg.TimeVec1024; cfg.TimeVec1024], BL_wind')),[],2);
BL_mean = mean(Data(:,:,BL_wind_idx(1):BL_wind_idx(2)),3);
if(strcmp(BL_corr,'covBL'))
    for ch = 1:128
        na_idx = find(~sum(isnan(BL_mean(:,ch)),2));
        X = BL_mean(na_idx,ch); %[ones(size(BL_mean,1),1) BL_mean(:,ch)];
        Data(na_idx,ch,:) = squeeze(Data(na_idx,ch,:)) - X*((X'*X)\X'*squeeze(Data(na_idx,ch,:)));
    end
elseif(strcmp(BL_corr,'minBL'))
    Data = bsxfun(@minus, Data, BL_mean);
end

% Fill in NaN Trials
for i = 1:length(Nan_trial_idx)
    Data = [Data(1:(Nan_trial_idx(i)-1),:,:); nan(1,size(Data,2),size(Data,3)); Data(Nan_trial_idx(i):end,:,:)];
end

% Sort Data
Data = Data(sort_idx,:,:);

% Save
Data_EEG = cfg;
Data_EEG.Data = Data;
Data_EEG.dim = 'trl_chan_time';

end