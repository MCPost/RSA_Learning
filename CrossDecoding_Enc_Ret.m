%% CrossDecoding Encoding and Retrieval


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs

% Subject Names
load('RSA_Data_Enc','Subj_names')

sub = 1;

load(['Preproc_EEG_Data/Encoding_object_locked/Enc_EEG_Data_',Subj_names{sub}])
eeg_cfg = [];
eeg_cfg.Name = Subj_names{sub};
eeg_cfg.chan_label = elecs;
eeg_cfg.chan_idx = ROI_all_idx';
eeg_cfg.Art_corr = false;
eeg_cfg.BL_corr = 'minBL';
eeg_cfg.BL_wind = [-0.26 -0.05];
[preprocEEG_Enc] = preproc_EEGData(eeg_cfg, Enc_Data_EEG);


load(['Preproc_EEG_Data/Retrieval_response_locked/Ret_EEG_Data_',Subj_names{sub}])
eeg_cfg = [];
eeg_cfg.Name = Subj_names{sub};
eeg_cfg.chan_label = elecs;
eeg_cfg.chan_idx = ROI_all_idx';
eeg_cfg.Art_corr = false;
eeg_cfg.BL_corr = 'demean';
eeg_cfg.BL_wind = [];
[preprocEEG_Ret] = preproc_EEGData(eeg_cfg, Ret_Data_EEG);



