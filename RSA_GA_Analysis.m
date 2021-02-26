%% RSA Learning - All Subjects Analysis

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
save('Data_EEG', 'Subj_names')
save('RSA_Data', 'Subj_names', 'elecs', 'ROI_all_idx', 'ROI_occ', 'ROI_occ_idx', 'ROI_temp', 'ROI_temp_idx')

%% Loop over Participants to create RSA Matrices
RSA_Mat_OCC = cell(1,2); RSA_Mat_TMP = cell(1,2);
MDS_Mat_OCC = cell(1,2); MDS_Mat_TMP = cell(1,2);
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
    rsa_cfg.meas16 = 'LDA';
    rsa_cfg.MNN = true;
    rsa_cfg.Cktl_blank_rm = true;
    rsa_cfg.only16 = true;
    
    rsa_cfg.curROI = ROI_occ_idx;
    rsa_cfg.curROI_name = 'occipital';
    RSA_Data_OCC = createRSA(rsa_cfg);
    
    assignin('base',['RSA_Data_OCC_',Subj_names{sub}],RSA_Data_OCC); 
    save('RSA_Data',['RSA_Data_OCC_',Subj_names{sub}],'-append')
    clear(['RSA_Mat_OCC_',Subj_names{sub}])
    
    rsa_cfg.curROI = ROI_temp_idx;
    rsa_cfg.curROI_name = 'temporal';
    RSA_Data_TMP = createRSA(rsa_cfg);
    
    assignin('base',['RSA_Data_TMP_',Subj_names{sub}],RSA_Mat_TMP); 
    save('RSA_Data',['RSA_Data_TMP_',Subj_names{sub}],'-append')
    clear(['RSA_Mat_TMP_',Subj_names{sub}])
    
end



TimeVec = RSA_Data_OCC.TimeVec;
TrialInfo = RSA_Data_OCC.TrialInfo;

figure
time_plot = -0.2:0.1:1.5;
for sbp = 1:18
    [~,t_idx] = min(abs(TimeVec - time_plot(sbp)));
    subplot(3,6,sbp)
    imagesc(squeeze(RSA_Data_OCC.RSA_16{1}(:,:,t_idx)))
    colorbar; title(sprintf('%4.0f ms', time_plot(sbp)*1000))
    axis square
end
    
figure
[~,t_idx] = min(abs(TimeVec - 0.3));
imagesc(squeeze(RSA_Data_OCC.RSA_16{1}(:,:,t_idx))); colorbar
title(sprintf('RSA Matrix at %4.0f ms', time_plot(sbp)*1000))
axis square

    
    
%% Hypotheses Matrix
Perceptual_Mat = zeros(size(TrialInfo,1));
Semantic_Mat = zeros(size(TrialInfo,1));
for i = 1:size(TrialInfo,1)-1
    for j = (i+1):size(TrialInfo,1)
        if(j ~= size(TrialInfo,1) - (i - 1))
            if(TrialInfo{i,5} == 1 && TrialInfo{j,5} == 1)
                Perceptual_Mat(i,j) = 1;
            elseif(TrialInfo{i,5} == 2 && TrialInfo{j,5} == 2)
                Perceptual_Mat(i,j) = 2;
            else
                Perceptual_Mat(i,j) = -1;
            end

            if(TrialInfo{i,7} == 1 && TrialInfo{j,7} == 1)
                Semantic_Mat(i,j) = 1;
            elseif(TrialInfo{i,7} == 2 && TrialInfo{j,7} == 2)
                Semantic_Mat(i,j) = 2;
            else
                Semantic_Mat(i,j) = -1;
            end
        end
    end
end

figure
subplot(1,2,1)
imagesc(Perceptual_Mat); title('Perceptual Hypothesis Matrix')
axis square
subplot(1,2,2)
imagesc(Semantic_Mat); title('Semantic Hypothesis Matrix')
axis square

Perceptual_Mat = zeros(16);
Semantic_Mat = zeros(16);
for i = 1:16-1
    for j = (i+1):16
        if(j ~= 16 - (i - 1))
            if(TrialInfo{8*i,5} == 1 && TrialInfo{8*j,5} == 1)
                Perceptual_Mat(i,j) = 1;
            elseif(TrialInfo{8*i,5} == 2 && TrialInfo{8*j,5} == 2)
                Perceptual_Mat(i,j) = 2;
            else
                Perceptual_Mat(i,j) = -1;
            end

            if(TrialInfo{8*i,7} == 1 && TrialInfo{8*j,7} == 1)
                Semantic_Mat(i,j) = 1;
            elseif(TrialInfo{8*i,7} == 2 && TrialInfo{8*j,7} == 2)
                Semantic_Mat(i,j) = 2;
            else
                Semantic_Mat(i,j) = -1;
            end
        end
    end
end

figure
subplot(1,2,1)
imagesc(Perceptual_Mat); title('Perceptual Hypothesis Matrix')
axis square
subplot(1,2,2)
imagesc(Semantic_Mat); title('Semantic Hypothesis Matrix')
axis square

    
    
%% Create RSA Time Courses

load RSA_Data rsa_cfg_tmp_CM
cfg = rsa_cfg_tmp_CM;
TimeVec = cfg.TimeVec;
TrialInfo = cfg.TrialInfo;

n = 2;
m = 5;

for sub = 1:length(Subj_names)

    %load('RSA_Data',['RSA_Mat_OCC_',Subj_names{sub}], ['MDS_Mat_OCC_CM'], ...
    %    ['RSA_Mat_TMP_',Subj_names{sub}], ['MDS_Mat_TMP_',Subj_names{sub}]);

    a = load('RSA_Data',['RSA_Mat_OCC_',Subj_names{sub}]);
    RSA_Mat_OCC = a.(['RSA_Mat_OCC_',Subj_names{sub}]);
    a = load('RSA_Data',['RSA_Mat_TMP_',Subj_names{sub}]);
    RSA_Mat_TMP = a.(['RSA_Mat_TMP_',Subj_names{sub}]);

    % Perceptual Dimension
    cur_data = zeros(size(RSA_Mat_OCC{n}{m},1));
    for tp = 1:length(TimeVec)
        cur_data = RSA_Mat_OCC{n}{m}(:,:,tp);
        RSA_Time.Perceptual_drawing(sub,tp) = nanmean(cur_data(Perceptual_Mat == 1));
        RSA_Time.Perceptual_picture(sub,tp) = nanmean(cur_data(Perceptual_Mat == 2));
        RSA_Time.Perceptual_within(sub,tp) = nanmean(cur_data(Perceptual_Mat > 0));
        RSA_Time.Perceptual_between(sub,tp) = nanmean(cur_data(Perceptual_Mat < 0));
    end

    % Semantic Dimension
    cur_data = zeros(size(RSA_Mat_OCC{n}{m},1));
    for tp = 1:length(TimeVec)
        cur_data = RSA_Mat_OCC{n}{m}(:,:,tp);
        RSA_Time.Semantic_animate(sub,tp) = nanmean(cur_data(Semantic_Mat == 1));
        RSA_Time.Semantic_inanimate(sub,tp) = nanmean(cur_data(Semantic_Mat == 2));
        RSA_Time.Semantic_within(sub,tp) = nanmean(cur_data(Semantic_Mat > 0));
        RSA_Time.Semantic_between(sub,tp) = nanmean(cur_data(Semantic_Mat < 0));
    end


    % Perceptual Dimension
    cur_data = zeros(size(RSA_Mat_TMP{n}{m},1));
    for tp = 1:length(TimeVec)
        cur_data = RSA_Mat_OCC{n}{m}(:,:,tp);
        RSA_Time.Perceptual_drawing(sub,tp) = nanmean(cur_data(Perceptual_Mat == 1));
        RSA_Time.Perceptual_picture(sub,tp) = nanmean(cur_data(Perceptual_Mat == 2));
        RSA_Time.Perceptual_within(sub,tp) = nanmean(cur_data(Perceptual_Mat > 0));
        RSA_Time.Perceptual_between(sub,tp) = nanmean(cur_data(Perceptual_Mat < 0));
    end

    % Semantic Dimension
    cur_data = zeros(size(RSA_Mat_TMP{n}{m},1));
    for tp = 1:length(TimeVec)
        cur_data = RSA_Mat_OCC{n}{m}(:,:,tp);
        RSA_Time.Semantic_animate(sub,tp) = nanmean(cur_data(Semantic_Mat == 1));
        RSA_Time.Semantic_inanimate(sub,tp) = nanmean(cur_data(Semantic_Mat == 2));
        RSA_Time.Semantic_within(sub,tp) = nanmean(cur_data(Semantic_Mat > 0));
        RSA_Time.Semantic_between(sub,tp) = nanmean(cur_data(Semantic_Mat < 0));
    end

end   



figure
subplot(2,2,1)
plot(TimeVec, nanmean(RSA_Time.Perceptual_drawing,1))
hold on
plot(TimeVec, nanmean(RSA_Time.Perceptual_picture,1))
hold off
ylabel('Corr'); xlabel('Time (s)'); title('Perceptual Dim Drawing vs. Picture')
xlim([-0.2 1.5]); legend('drawing','picture')
subplot(2,2,2)
plot(TimeVec, nanmean(RSA_Time.Perceptual_within,1))
hold on
plot(TimeVec, nanmean(RSA_Time.Perceptual_between,1))
hold off
ylabel('Corr'); xlabel('Time (s)'); title('Perceptual Dim Within vs. Between Cat')
xlim([-0.2 1.5]); legend('perceptual within','perceptual between')
subplot(2,2,3)
plot(TimeVec, nanmean(RSA_Time.Semantic_animate,1))
hold on
plot(TimeVec, nanmean(RSA_Time.Semantic_inanimate,1))
hold off
ylabel('Corr'); xlabel('Time (s)'); title('Semantic Dim Animate vs. Inanimate')
xlim([-0.2 1.5]); legend('Animate','Inanimate')
subplot(2,2,4)
plot(TimeVec, nanmean(RSA_Time.Semantic_within,1))
hold on
plot(TimeVec, nanmean(RSA_Time.Semantic_between,1))
hold off
ylabel('Corr'); xlabel('Time (s)'); title('Semantic Dim Within vs. Between Cat')
xlim([-0.2 1.5]); legend('semantic within','semantic between')

figure
plot(TimeVec, nanmean(RSA_Time.Perceptual_within,1))
hold on
plot(TimeVec, nanmean(RSA_Time.Semantic_within,1))
hold off
ylabel('Corr'); xlabel('Time (s)'); title('Perceptual vs. Semantic')
xlim([-0.2 1.5]); legend('Perceptual','Semantic')



