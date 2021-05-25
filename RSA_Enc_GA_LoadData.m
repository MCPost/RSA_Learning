%% Load in RSA Encoding Grand Average Data

function RSA_Enc_GA_LoadData(sub, Subj_names, ROI_idx, save_name)


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs


%% Initialize Mat files
%save('RSA_Data_Enc', 'Subj_names', 'elecs', 'ROI_all_idx')

%% Loop over Participants to create RSA Matrices
%for sub = 1:length(Subj_names)
    
    %% Import Data
    load(['Preproc_EEG_Data/Encoding_object_locked/Enc_EEG_Data_',Subj_names{sub}])
    eeg_cfg = [];
    eeg_cfg.Name = Subj_names{sub};
    eeg_cfg.chan_label = elecs;
    eeg_cfg.chan_idx = ROI_all_idx'; 
    eeg_cfg.Art_corr = false;
    eeg_cfg.BL_corr = 'minBL';
    eeg_cfg.BL_wind = [-0.26 -0.05];
    [Data_EEG] = preproc_EEGData(eeg_cfg, Enc_Data_EEG);
    
    %assignin('base',['Data_EEG_',Subj_names{sub}],Data_EEG)
    %save('Data_EEG',['Data_EEG_',Subj_names{sub}],'-append')
    %clear(['Data_EEG_',Subj_names{sub}])
    
    
    %% RSA Function
    rsa_cfg = Data_EEG;
    rsa_cfg.slide_window = 0.080; 
    rsa_cfg.slide_step = 0.004;
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
