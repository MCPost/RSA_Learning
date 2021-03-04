%% Load in RSA Grand Average Data


%% Subject Names
Subj_names = {'AG','CEWD','CM','DS','FVM'};


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
elecs = {'A1';'A2';'A3';'A4';'A5';'A6';'A7';'A8';'A9';'A10';'A11';'A12';'A13';'A14';'A15';'A16';'A17';...
          'A18';'A19';'A20';'A21';'A22';'A23';'A24';'A25';'A26';'A27';'A28';'A29';'A30';'A31';'A32';'B1';...
          'B2';'B3';'B4';'B5';'B6';'B7';'B8';'B9';'B10';'B11';'B12';'B13';'B14';'B15';'B16';'B17';'B18';...
          'B19';'B20';'B21';'B22';'B23';'B24';'B25';'B26';'B27';'B28';'B29';'B30';'B31';'B32';'C1';'C2';...
          'C3';'C4';'C5';'C6';'C7';'C8';'C9';'C10';'C11';'C12';'C13';'C14';'C15';'C16';'C17';'C18';'C19';...
          'C20';'C21';'C22';'C23';'C24';'C25';'C26';'C27';'C28';'C29';'C30';'C31';'C32';'D1';'D2';'D3';...
          'D4';'D5';'D6';'D7';'D8';'D9';'D10';'D11';'D12';'D13';'D14';'D15';'D16';'D17';'D18';'D19';'D20';...
          'D21';'D22';'D23';'D24';'D25';'D26';'D27';'D28';'D29';'D30';'D31';'D32'};

% All Electrodes
ROI_all_idx = 1:128;

% Occipital Cortex ROI
ROI_occ = {'A8','A9','A10','A11','A12','A13','A14','A15','A16','A17','A21','A22','A23','A24','A25','A26',...
           'A27','A28','A29','A30','B5','B6','B7','B8','B9',...
           'D32','D31','D30','D29','A6','A7','A5','A18','A19','A20','A32','A31','B3','B4','B13','B12','B11','B10'};
ROI_occ_idx = find(cell2mat(cellfun(@(x) any(strcmp(x, ROI_occ)), elecs, 'UniformOutput', 0)));

% Temporal Cortex ROI
ROI_temp = {'B18','B17','B16','B15','B14','B22','B23','B24','B25','B26','B31','B30','B29','B28','B27','C3',...
            'C4','C5','C6','C7','D28','D27','D26','D25','D24','D19','D20','D21','D22','D23','D8','D9','D10',...
            'D11','D12','D7','D6','D5','D4','D3'};
ROI_temp_idx = find(cell2mat(cellfun(@(x) any(strcmp(x, ROI_temp)), elecs, 'UniformOutput', 0)));


%% Initialize Mat files
save('RSA_Data', 'Subj_names', 'elecs', 'ROI_all_idx')

%% Loop over Participants to create RSA Matrices
for sub = 1:length(Subj_names)
    
    %% Import Data
    eeg_cfg = [];
    eeg_cfg.Name = Subj_names{sub};
    eeg_cfg.chan_label = elecs;
    eeg_cfg.chan_idx = ROI_occ_idx'; 
    eeg_cfg.Art_corr = true;
    eeg_cfg.BL_corr = 'minBL';
    eeg_cfg.BL_wind = [-0.26 -0.05];
    [Data_EEG] = createEEGData(eeg_cfg);
    
    %assignin('base',['Data_EEG_',Subj_names{sub}],Data_EEG)
    %save('Data_EEG',['Data_EEG_',Subj_names{sub}],'-append')
    %clear(['Data_EEG_',Subj_names{sub}])
    
    
    %% RSA Function
    rsa_cfg = Data_EEG;
    rsa_cfg.slide_window = 0.080; 
    rsa_cfg.slide_step = 0.004;
    rsa_cfg.window_average = 'gaussian';
    rsa_cfg.meas128 = '';
    rsa_cfg.meas16 = 'SVM';
    rsa_cfg.MNN = true;
    rsa_cfg.Cktl_blank_rm = true;
    rsa_cfg.only16 = true;
    
    rsa_cfg.curROI = ROI_occ_idx;
    rsa_cfg.curROI_name = 'occipital';
    RSA_Data_OCC.(['RSA_Data_OCC_',Subj_names{sub}]) = createRSA(rsa_cfg);
    
    save('RSA_Data_SVM','-struct','RSA_Data_OCC','-append')
    clear('RSA_Mat_OCC')
    
    rsa_cfg.curROI = ROI_temp_idx;
    rsa_cfg.curROI_name = 'temporal';
    RSA_Data_TMP.(['RSA_Data_TMP_',Subj_names{sub}]) = createRSA(rsa_cfg);
    
    save('RSA_Data_SVM','-struct','RSA_Data_TMP','-append')
    clear('RSA_Mat_TMP')
    
end


