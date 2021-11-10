%% Load in RSA Retrieval Grand Average Data

function RSA_Ret_GA_LoadData(sub, Subj_names, ROI_idx, save_name)


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs


%% Initialize Mat files
%save('RSA_Data_Ret', 'Subj_names', 'elecs', 'ROI_all_idx')

%% Loop over Participants to create RSA Matrices
%for sub = 1:length(Subj_names)
    
    %% Import Data
    load(['Preproc_EEG_Data/Retrieval_response_locked/Ret_EEG_Data_',Subj_names{sub}])
    eeg_cfg = [];
    eeg_cfg.Name = Subj_names{sub};
    eeg_cfg.chan_label = elecs;
    eeg_cfg.chan_idx = ROI_all_idx; %ROI_occ_idx'; 
    eeg_cfg.Art_corr = false;
    eeg_cfg.BL_corr = 'demean';
    eeg_cfg.BL_wind = [-4 0.2];
    [Data_EEG] = preproc_EEGData(eeg_cfg, Ret_Data_EEG);
    
    
    %% RSA Function
    rsa_cfg = Data_EEG;
    rsa_cfg.slide_window = 0.080; 
    rsa_cfg.slide_step = 0.008;
    rsa_cfg.fwhm = 0.040;
    rsa_cfg.window_average = 'gaussian';
    rsa_cfg.meas128 = '';
    rsa_cfg.meas16 = {'LDA', 'SVM', 'euclidian'};
    rsa_cfg.MNN = true;
    rsa_cfg.Cktl_blank_rm = false;
    rsa_cfg.only16 = true;
    rsa_cfg.noMDS = true;
    rsa_cfg.labelshuffperm = false;
    rsa_cfg.n_perms = 1000;
    
    for r = 1:size(ROI_idx,1)
        rsa_cfg.curROI = ROI_idx{r,1};
        rsa_cfg.curROI_name = ROI_idx{r,2};
        RSA_Data.(['RSA_Data_',Subj_names{sub}]).(ROI_idx{r,1}) = createRSA(rsa_cfg);
    end
    
    save(save_name,'-struct','RSA_Data','-append')
    clear('RSA_Data')
    
%end

end


%Data = Ret_Data_EEG.Retrieval_Data(cell2mat(Ret_Data_EEG.TrialInfo(cell2mat(Ret_Data_EEG.TrialInfo(:,1)) > 0,2)),:);

%hist(cell2mat(Data(:,9)),20)





