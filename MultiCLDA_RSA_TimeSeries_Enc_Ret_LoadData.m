%% LDA RSA Time Course Encoding and Retrieval Load Data


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs 
ROI = {'OCC';'TMP';'FRT';'PRT'};
ROI_idx = {ROI_occ_idx; ROI_temp_idx; ROI_front_idx; ROI_pari_idx};

% Hypothesis Matrix
Hypothesis_Matrix

% Subject Names
load('RSA_Data_Enc','Subj_names')

% Gaussian Average across time window and downsample
window_average  = 'gaussian';
slide_window_s  = 0.080; %s
slide_step_s    = 0.010; %s
fwhm            = 0.040; %s
if(strcmp(window_average,'gaussian'))
    average_kern = @(x, dim, gsize, sd) sum(bsxfun(@times, x, (gausswin(gsize,sd)./sum(gausswin(gsize,sd)))'),dim);
    av = 'gaussian'; % Gaussian Average
else
    average_kern = @(x, dim, gsize, sd) mean(x,dim);
    av = 'uniform';  % Uniform Average
end

samplingrate = 1024; %Hz

% Matrix Indices to choose
UTmat = triu(ones(16)).*~eye(16);
UTM_idx = find(UTmat(:) == 1);


MultiCLDA_RSA_Enc_Data = [];
MultiCLDA_RSA_Enc_Data.window_average = window_average;
MultiCLDA_RSA_Enc_Data.slide_window = slide_window_s;
MultiCLDA_RSA_Enc_Data.slide_step = slide_step_s;
MultiCLDA_RSA_Enc_Data.fwhm = fwhm;
MultiCLDA_RSA_Enc_Data.ROI = ROI;
MultiCLDA_RSA_Enc_Data.TimeEnc = [];
[~,MultiCLDA_RSA_Enc_Data.Cond_idx.Perc_WI_BT_idx(:,1),~] = intersect(UTM_idx,find(Perceptual_Mat_red16(:) > 0));
[~,MultiCLDA_RSA_Enc_Data.Cond_idx.Perc_WI_BT_idx(:,2),~] = intersect(UTM_idx,find(Perceptual_Mat_red16(:) < 0));
[~,MultiCLDA_RSA_Enc_Data.Cond_idx.Sem_WI_BT_idx(:,1),~]  = intersect(UTM_idx,find(Semantic_Mat_red16(:) > 0));
[~,MultiCLDA_RSA_Enc_Data.Cond_idx.Sem_WI_BT_idx(:,2),~]  = intersect(UTM_idx,find(Semantic_Mat_red16(:) < 0));


MultiCLDA_RSA_Ret_Data = [];
MultiCLDA_RSA_Ret_Data.window_average = window_average;
MultiCLDA_RSA_Ret_Data.slide_window = slide_window_s;
MultiCLDA_RSA_Ret_Data.slide_step = slide_step_s;
MultiCLDA_RSA_Ret_Data.fwhm = fwhm;
MultiCLDA_RSA_Ret_Data.ROI = ROI;
MultiCLDA_RSA_Ret_Data.TimeRet = [];
[~,MultiCLDA_RSA_Ret_Data.Cond_idx.Perc_WI_BT_idx(:,1),~] = intersect(UTM_idx,find(Perceptual_Mat_red16(:) > 0));
[~,MultiCLDA_RSA_Ret_Data.Cond_idx.Perc_WI_BT_idx(:,2),~] = intersect(UTM_idx,find(Perceptual_Mat_red16(:) < 0));
[~,MultiCLDA_RSA_Ret_Data.Cond_idx.Sem_WI_BT_idx(:,1),~]  = intersect(UTM_idx,find(Semantic_Mat_red16(:) > 0));
[~,MultiCLDA_RSA_Ret_Data.Cond_idx.Sem_WI_BT_idx(:,2),~]  = intersect(UTM_idx,find(Semantic_Mat_red16(:) < 0));


fprintf('\n')
nbytes = fprintf('Subject: %s  --  Progress Cross Decoding:  0.0 %%',Subj_names{sub});
for sub = 1:length(Subj_names)

    % Load Encoding Data
    load(['Preproc_EEG_Data/Encoding_object_locked/Enc_EEG_Data_',Subj_names{sub}])
    eeg_cfg = [];
    eeg_cfg.Name = Subj_names{sub};
    eeg_cfg.chan_label = elecs;
    eeg_cfg.chan_idx = ROI_all_idx';
    eeg_cfg.Art_corr = false;
    eeg_cfg.BL_corr = 'minBL';
    eeg_cfg.BL_wind = [-0.26 -0.05];
    [preprocEEG_Enc] = preproc_EEGData(eeg_cfg, Enc_Data_EEG);

    % Load Retrieval Data
    load(['Preproc_EEG_Data/Retrieval_response_locked/Ret_EEG_Data_',Subj_names{sub}])
    eeg_cfg = [];
    eeg_cfg.Name = Subj_names{sub};
    eeg_cfg.chan_label = elecs;
    eeg_cfg.chan_idx = ROI_all_idx';
    eeg_cfg.Art_corr = false;
    eeg_cfg.BL_corr = 'demean';
    eeg_cfg.BL_wind = [-4 0.2];
    [preprocEEG_Ret] = preproc_EEGData(eeg_cfg, Ret_Data_EEG);

    % Downsampling
    slide_window_idx = round(samplingrate*(slide_window_s/2));
    slide_step_idx = round(samplingrate*slide_step_s);

    TimeVec_idx = 1:slide_step_idx:length(preprocEEG_Enc.TimeVec1024);
    TimeVecEnc = preprocEEG_Enc.TimeVec1024(TimeVec_idx);
    preprocEEG_Enc.Data_dwnsmp = zeros(size(preprocEEG_Enc.Data,1), size(preprocEEG_Enc.Data,2), length(TimeVecEnc));
    for tp = 1:length(TimeVecEnc) 
        time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, length(preprocEEG_Enc.TimeVec1024));
        preprocEEG_Enc.Data_dwnsmp(:,:,tp) = squeeze(average_kern(permute(preprocEEG_Enc.Data(:,:,time_window),[1 3 2]), 2, length(time_window), length(time_window)/round(samplingrate*(fwhm/2))*.44));
    end

    TimeVec_idx = 1:slide_step_idx:length(preprocEEG_Ret.TimeVec1024);
    TimeVecRet = preprocEEG_Ret.TimeVec1024(TimeVec_idx);
    preprocEEG_Ret.Data_dwnsmp = zeros(size(preprocEEG_Ret.Data,1), size(preprocEEG_Ret.Data,2), length(TimeVecRet));
    for tp = 1:length(TimeVecRet) 
        time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, length(preprocEEG_Ret.TimeVec1024));
        preprocEEG_Ret.Data_dwnsmp(:,:,tp) = squeeze(average_kern(permute(preprocEEG_Ret.Data(:,:,time_window),[1 3 2]), 2, length(time_window), length(time_window)/round(samplingrate*(fwhm/2))*.44));
    end

    % Save TimeVecs
    MultiCLDA_RSA_Enc_Data.TimeEnc = TimeVecEnc;
    MultiCLDA_RSA_Ret_Data.TimeRet = TimeVecRet;
    
    for r = 1:size(ROI,1)
    
        % Initialize Data
        if(sub == 1)
            MultiCLDA_RSA_Enc_Data.(ROI{r}).Data = zeros(length(Subj_names),length(TimeVecEnc),length(UTM_idx));
            MultiCLDA_RSA_Ret_Data.(ROI{r}).Data = zeros(length(Subj_names),length(TimeVecRet),length(UTM_idx));
            MultiCLDA_RSA_Enc_Data.(ROI{r}).Sensor_weights = zeros(length(Subj_names),length(TimeVecEnc),length(ROI_idx{r}));
            MultiCLDA_RSA_Ret_Data.(ROI{r}).Sensor_weights = zeros(length(Subj_names),length(TimeVecRet),length(ROI_idx{r}));
        end


        %% Crossdecoding

        % Current Data
        CurData_Enc = preprocEEG_Enc.Data_dwnsmp(~isnan(sum(preprocEEG_Enc.Data_dwnsmp(:,:,1),2)),ROI_idx{r},:);
        CurData_Ret = preprocEEG_Ret.Data_dwnsmp(~isnan(sum(preprocEEG_Ret.Data_dwnsmp(:,:,1),2)),ROI_idx{r},:);

        % Perceptual Dimension
        Cond_Enc = cell2mat(preprocEEG_Enc.TrialInfo(cell2mat(preprocEEG_Enc.TrialInfo(:,1)) > 0,9)) + (cell2mat(preprocEEG_Enc.TrialInfo(cell2mat(preprocEEG_Enc.TrialInfo(:,1)) > 0,5)) == 2).*8;
        Cond_Ret = cell2mat(preprocEEG_Ret.TrialInfo(cell2mat(preprocEEG_Ret.TrialInfo(:,1)) > 0,9)) + (cell2mat(preprocEEG_Ret.TrialInfo(cell2mat(preprocEEG_Ret.TrialInfo(:,1)) > 0,5)) == 2).*8;

        cfg = [];
        cfg.metric      = 'confusion';
        cfg.classifier  = 'multiclass_lda';
        %cfg.cv          = 'holdout';
        %cfg.p           = 0.3;
        cfg.cv          = 'kfold';
        cfg.k           = 4;
        cfg.repeat      = 5;
        cfg.feedback    = 1;
        
        [perf,~,~,weights] = mv_classify_weights(cfg, CurData_Enc, Cond_Enc);
        for tp = 1:size(perf,1)
            dgn_vec = diag(squeeze(perf(tp,:,:))) + double(diag(squeeze(perf(tp,:,:))) == 0);
            S = sqrt(bsxfun(@rdivide, squeeze(perf(tp,:,:)), dgn_vec) * bsxfun(@rdivide, squeeze(perf(tp,:,:)), dgn_vec)');
            MultiCLDA_RSA_Enc_Data.(ROI{r}).Data(sub,tp,:) = 1 - S(UTM_idx)';
            MultiCLDA_RSA_Enc_Data.(ROI{r}).Sensor_weights(sub,tp,:) = (cov(CurData_Enc(:,:,tp))*mean(mean(weights(:,:,:,tp),3),2))' / cov(CurData_Enc(:,:,tp)*mean(mean(weights(:,:,:,tp),3),2));
        end
        
        [perf,~,~,weights] = mv_classify_weights(cfg, CurData_Ret, Cond_Ret);
        for tp = 1:size(perf,1)
            dgn_vec = diag(squeeze(perf(tp,:,:))) + double(diag(squeeze(perf(tp,:,:))) == 0);
            S = sqrt(bsxfun(@rdivide, squeeze(perf(tp,:,:)), dgn_vec) * bsxfun(@rdivide, squeeze(perf(tp,:,:)), dgn_vec)');
            MultiCLDA_RSA_Ret_Data.(ROI{r}).Data(sub,tp,:) = 1 - S(UTM_idx)';
            MultiCLDA_RSA_Ret_Data.(ROI{r}).Sensor_weights(sub,tp,:) = (cov(CurData_Ret(:,:,tp))*mean(mean(weights(:,:,:,tp),3),2))' / cov(CurData_Ret(:,:,tp)*mean(mean(weights(:,:,:,tp),3),2));
        end
        
        
        fprintf(repmat('\b',1,nbytes))
        nbytes = fprintf('Subject: %s  --  Progress Cross Decoding:  %3.1f %%',Subj_names{sub},(r / size(ROI,1))*100);
    end
    fprintf('\n')
end
fprintf('\n')


%% Save Data

save('MultiCLDA_RSA_Enc','MultiCLDA_RSA_Ret_Enc')
save('MultiCLDA_RSA_Ret','MultiCLDA_RSA_Ret_Data')

