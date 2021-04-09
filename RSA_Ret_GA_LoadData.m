%% Load in RSA Retrieval Grand Average Data

function RSA_Ret_GA_LoadData(sub, Subj_names)


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
    eeg_cfg.BL_wind = [-0.26 -0.05];
    [Data_EEG] = preproc_EEGData(eeg_cfg, Ret_Data_EEG);
    
    
    %% RSA Function
    rsa_cfg = Data_EEG;
    rsa_cfg.slide_window = 0.080; 
    rsa_cfg.slide_step = 0.008;
    rsa_cfg.window_average = 'gaussian';
    rsa_cfg.meas128 = '';
    rsa_cfg.meas16 = {'LDA', 'SVM', 'euclidian', 'euclidian w.c.c.'};
    rsa_cfg.MNN = true;
    rsa_cfg.Cktl_blank_rm = true;
    rsa_cfg.only16 = true;
    rsa_cfg.noMDS = true;
    
    rsa_cfg.curROI = ROI_occ_idx;
    rsa_cfg.curROI_name = 'occipital';
    RSA_Data.(['RSA_Data_',Subj_names{sub}]).OCC = createRSA(rsa_cfg);
    
    rsa_cfg.curROI = ROI_temp_idx;
    rsa_cfg.curROI_name = 'temporal';
    RSA_Data.(['RSA_Data_',Subj_names{sub}]).TMP = createRSA(rsa_cfg);
    
    rsa_cfg.curROI = ROI_front_idx;
    rsa_cfg.curROI_name = 'frontal';
    RSA_Data.(['RSA_Data_',Subj_names{sub}]).FRT = createRSA(rsa_cfg);
    
    rsa_cfg.curROI = ROI_cent_idx;
    rsa_cfg.curROI_name = 'central';
    RSA_Data.(['RSA_Data_',Subj_names{sub}]).CNT = createRSA(rsa_cfg);
    
    rsa_cfg.curROI = ROI_pari_idx;
    rsa_cfg.curROI_name = 'parietal';
    RSA_Data.(['RSA_Data_',Subj_names{sub}]).PRT = createRSA(rsa_cfg);
    
    save('RSA_Data_Ret','-struct','RSA_Data','-append')
    clear('RSA_Data')
    
%end

end
