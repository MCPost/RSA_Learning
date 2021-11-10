%% CrossDecoding Encoding and Retrieval Load Data

function CrossDecoding_Enc_Ret_LoadData(msr, cfg, pair, save_name)

% Measure
cfg.msr = msr;

% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs 
cfg.ROI = ROI;
cfg.ROI_idx = {ROI_occ_idx; ROI_temp_idx; ROI_front_idx; ROI_pari_idx};

% Hypotheses Matrix
Hypothesis_Matrix


%% Matrix Indices to choose
UTmat = triu(ones(16)).*~eye(16);
UTM_idx = find(UTmat(:) == 1);
SensUTmat = UTmat(1:end-1,:);
SensUTM_idx = find(SensUTmat(:) == 1);


%% Load Data 

% Subject Names
load('RSA_Data_Enc','Subj_names')

for sub = 1:length(Subj_names)
    
    cfg.subs = Subj_names{sub};
    
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

    
    %% Crossdecode time courses for each subject

    if(strcmp(pair,'Enc - Enc'))
        Data1 = preprocEEG_Enc;
        Data2 = preprocEEG_Enc;
    elseif(strcmp(pair,'Ret - Ret'))
        Data1 = preprocEEG_Ret;
        Data2 = preprocEEG_Ret;
    else
        Data1 = preprocEEG_Enc; 
        Data2 = preprocEEG_Ret;
    end

    CrossDecoding.(['CrossDecoding_msr',num2str(msr),'_',Subj_names{sub}]) = create_xdecode(cfg, Data1, Data2);

    save([save_name,num2str(msr),'_',Subj_names{sub}],'-struct','CrossDecoding')

end


