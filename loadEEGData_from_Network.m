%% Load EEG Data from the Network

% Get Path to EEG data

eeg_fullpath = 'Z:/Christopher/RSA_for_Christopher/EGG_data_preprocessed/';
eeg_encpath = 'Full_clean_encoding_2/';
eeg_retpath = 'Full_clean_retrieval_RL_4seconds/';

if(exist(eeg_fullpath,'dir') ~= 7)
    error('Path does not exist. Network probably not mounted.')
end

% Get logfiles and Names
[~,message,~] = fileattrib('Logs\*');
currentdir = pwd;
filenames = strrep({message([message.directory] == 0).Name}',[currentdir,'\'],'');
filenames = regexprep(filenames,'\','/');

tok_names = unique(cellfun(@(x) x{1}, regexp(filenames,'_(\w*).','tokens','once'),'UniformOutput', 0));
Subj_names = tok_names(cellfun(@isempty, regexp(tok_names,'_')));

% leave out Subj 'TVV' 
Subj_names(strcmp(Subj_names,'TVV')) = [];


for sub = 1:length(Subj_names)
    
    Name = Subj_names{sub};
    
    
    %% Load Behavioral Data (Logfile)
    
    fprintf('\n')
    fprintf('Load Logfile Data for Subject %s!', Name)
    fprintf('\n')
    
    if(sub == 20)
        logfile = importdata(['Logs/resultfile_',Name,'_nocutted.csv']);
    else
        logfile = importdata(['Logs/resultfile_',Name,'.csv']);
    end
    
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

    
    %% Load EEG Data
    
    fprintf('Load Encoding EEG Data for Subject %s!', Name)
    fprintf('\n')
    
    Preproc_Data_1 = load([eeg_fullpath,eeg_encpath,'Prepro_2_encoding_full_',Name,'.mat'], ['r_1_',Name]);
    Preproc_Data_1 = Preproc_Data_1.(cell2mat(fieldnames(Preproc_Data_1)));
    Preproc_Data_2 = load([eeg_fullpath,eeg_encpath,'Prepro_2_encoding_full_',Name,'.mat'], ['r_2_',Name]);
    Preproc_Data_2 = Preproc_Data_2.(cell2mat(fieldnames(Preproc_Data_2)));
    
    
    %% Marry Behavioral and EEG Data
    
    fprintf('Marry Encoding Behavioral and EEG Data for Subject %s!', Name)
    fprintf('\n')
    
    TrialInfo = cell(size(Encoding_Data,1),12);
    TI2temp = num2cell(Preproc_Data_2.trialinfo);
    draw_idx = find(cell2mat(Encoding_Data(:,6)) == 1); cnt = 1;
    for i = 1:length(draw_idx)
        if(ismember(draw_idx(i), cell2mat(TI2temp(:,2))))
            TrialInfo(i,:) = [TI2temp(cnt,1:4) Encoding_Data(draw_idx(i),[6 14 16 15 7 17 5 4])];
            cnt = cnt + 1;
        else
            TrialInfo(i,:) = [{0}, {draw_idx(i)} Encoding_Data(draw_idx(i),[6 7 6 14 16 15 7 17 5 4])];
        end
    end
    TI1temp = num2cell(Preproc_Data_1.trialinfo);
    pic_idx = find(cell2mat(Encoding_Data(:,6)) == 2); cnt = 1;
    for i = 1:length(pic_idx)
        if(ismember(pic_idx(i), cell2mat(TI1temp(:,2))))
            TrialInfo(64+i,:) = [TI1temp(cnt,1:4) Encoding_Data(pic_idx(i),[6 14 16 15 7 17 5 4])];
            cnt = cnt + 1;
        else
            TrialInfo(64+i,:) = [{0}, {pic_idx(i)} Encoding_Data(pic_idx(i),[6 7 6 14 16 15 7 17 5 4])];
        end
    end

    % Get Index of NaN trials
    Nan_trial_idx = find(cell2mat(TrialInfo(:,1)) == 0);


    %% Create EEG Data Struct for Encoding
    [~,sort_idx] = sortrows(cell2mat(TrialInfo(:,[5 7 9])));

    cfg.TrialInfo = TrialInfo(sort_idx,:);
    cfg.TimeVec1024 = Preproc_Data_1.time;

    Data = [Preproc_Data_2.trial; Preproc_Data_1.trial];
    
    % Fill in NaN Trials
    for i = 1:length(Nan_trial_idx)
        Data = [Data(1:(Nan_trial_idx(i)-1),:,:); nan(1,size(Data,2),size(Data,3)); Data(Nan_trial_idx(i):end,:,:)];
    end

    % Sort Data
    Data = Data(sort_idx,:,:);

    % Save
    Enc_Data_EEG = cfg;
    Enc_Data_EEG.Data = Data;
    Enc_Data_EEG.dim = 'trl_chan_time';
    
    fprintf('Save Encoding Data for Subject %s!', Name)
    fprintf('\n')
    save(['Preproc_EEG_Data/Encoding_object_locked/Enc_EEG_Data_',Name], 'Enc_Data_EEG')
    
    clear cfg Enc_Data_EEG Data Preproc_Data_1 Preproc_Data_2 Encoding_Data Nan_trial_idx sort_idx TrialInfo
    
    
    %% Get Retrieval from logfile and Categories of Objects

    names = {'block_number','trial_retrieval','trial_total_ret','object_name','object_number','object_percep',... 
             'object_category','cue_retrieval','RT_reinstatement','catch_trial_type','acc_catch','RT_catch_trial',...
             'catch_trial_2_type','acc_catch_2','RT_catch__2_trial','total_responses_distractor','accuracy_distractor',...
             'OnsetRetrievalBlockStart','CueRetrievalOnset','ResponseRetrievalOnset','CatchTrialOnset',...
             'ResponseCatchOnset','CatchTrialOnset_2','ResponseCatchOnset_2','episode_remembered'};

    Retrieval_Data = logfile(strcmp(logfile(:,2),'retrieval'), cell2mat(cellfun(@(x) ismember(x,names), strtrim(logfile(1,:)),'UniformOutput', 0)));
    % block_number, trial_encoding, trial_total_enc, object_name,
    % object_number, object_percep, object_category, RT_encoding,
    % cue_retrieval, OnsetEncodingBlockStart, WordOnsetEncoding, 
    % ObjectContextOnset, ResponseEncodingOnset, PerceptualDim, SemanticDim,
    % Category

    Retrieval_Data(cell2mat(Retrieval_Data(:,6)) == 1,26) = {'Drawing'};
    Retrieval_Data(cell2mat(Retrieval_Data(:,6)) == 2,26) = {'Picture'};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) < 5,27) = {'Animate'};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) > 4,27) = {'Inanimate'};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) < 5,28) = {1};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) > 4,28) = {2};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) == 1,29) = {'insect'};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) == 2,29) = {'bird'};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) == 3,29) = {'mammal'};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) == 4,29) = {'sea'};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) == 5,29) = {'fruit'};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) == 6,29) = {'electronic'};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) == 7,29) = {'veggie'};
    Retrieval_Data(cell2mat(Retrieval_Data(:,7)) == 8,29) = {'clothe'};
    cfg.Retrieval_Data = Retrieval_Data;

    
    %% Load EEG Data
    
    fprintf('Load Retrieval EEG Data for Subject %s!', Name)
    fprintf('\n')
    
    Preproc_Data_1 = load([eeg_fullpath,eeg_retpath,'Prepro_retrieval_RL_4seconds_full_',Name,'.mat'], ['r_1_',Name]);
    Preproc_Data_1 = Preproc_Data_1.(cell2mat(fieldnames(Preproc_Data_1)));
    Preproc_Data_2 = load([eeg_fullpath,eeg_retpath,'Prepro_retrieval_RL_4seconds_full_',Name,'.mat'], ['r_2_',Name]);
    Preproc_Data_2 = Preproc_Data_2.(cell2mat(fieldnames(Preproc_Data_2)));
    
    
    %% Marry Behavioral and EEG Data
    
    fprintf('Marry Retrieval Behavioral and EEG Data for Subject %s!', Name)
    fprintf('\n')
    
    TrialInfo = cell(size(Retrieval_Data,1),12);
    TI2temp = num2cell(Preproc_Data_2.trialinfo);
    draw_idx = find(cell2mat(Retrieval_Data(:,6)) == 1); cnt = 1;
    for i = 1:length(draw_idx)
        if(ismember(draw_idx(i), cell2mat(TI2temp(:,2))))
            TrialInfo(i,:) = [TI2temp(cnt,1:4) Retrieval_Data(draw_idx(i),[6 26 28 27 7 29 5 4])];
            cnt = cnt + 1;
        else
            TrialInfo(i,:) = [{0}, {draw_idx(i)} Retrieval_Data(draw_idx(i),[6 7 6 26 28 27 7 29 5 4])];
        end
    end
    TI1temp = num2cell(Preproc_Data_1.trialinfo);
    pic_idx = find(cell2mat(Retrieval_Data(:,6)) == 2); cnt = 1;
    for i = 1:length(pic_idx)
        if(ismember(pic_idx(i), cell2mat(TI1temp(:,2))))
            TrialInfo(64+i,:) = [TI1temp(cnt,1:4) Retrieval_Data(pic_idx(i),[6 26 28 27 7 29 5 4])];
            cnt = cnt + 1;
        else
            TrialInfo(64+i,:) = [{0}, {pic_idx(i)} Retrieval_Data(pic_idx(i),[6 7 6 26 28 27 7 29 5 4])];
        end
    end

    % Get Index of NaN trials
    Nan_trial_idx = find(cell2mat(TrialInfo(:,1)) == 0);


    %% Create EEG Data Struct for Retrieval
    [~,sort_idx] = sortrows(cell2mat(TrialInfo(:,[5 7 9])));

    cfg.TrialInfo = TrialInfo(sort_idx,:);
    cfg.TimeVec1024 = Preproc_Data_1.time;

    Data = [Preproc_Data_2.trial; Preproc_Data_1.trial];
    
    % Fill in NaN Trials
    for i = 1:length(Nan_trial_idx)
        Data = [Data(1:(Nan_trial_idx(i)-1),:,:); nan(1,size(Data,2),size(Data,3)); Data(Nan_trial_idx(i):end,:,:)];
    end

    % Sort Data
    Data = Data(sort_idx,:,:);

    % Save
    Ret_Data_EEG = cfg;
    Ret_Data_EEG.Data = Data;
    Ret_Data_EEG.dim = 'trl_chan_time';
    
    fprintf('Save Retrieval Data for Subject %s!', Name)
    fprintf('\n')
    save(['Preproc_EEG_Data/Retrieval_response_locked/Ret_EEG_Data_',Name], 'Ret_Data_EEG')
    
    clear cfg Ret_Data_EEG Data Preproc_Data_1 Preproc_Data_2 Retrieval_Data Nan_trial_idx sort_idx TrialInfo
    
end

