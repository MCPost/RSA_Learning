%% RSA Learning - Single Subject Analysis Example

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
           'A27','A28','A29','A30','B5','B6','B7','B8','B9'};
ROI_occ_idx = find(cell2mat(cellfun(@(x) any(strcmp(x, ROI_occ)), elecs, 'UniformOutput', 0)));

% Temporal Cortex ROI
ROI_temp = {'B18','B17','B16','B15','B14','B22','B23','B24','B25','B26','B31','B30','B29','B28','B27','C3',...
            'C4','C5','C6','C7','D28','D27','D26','D25','D24','D19','D20','D21','D22','D23','D8','D9','D10',...
            'D11','D12','D7','D6','D5','D4','D3'};
ROI_temp_idx = find(cell2mat(cellfun(@(x) any(strcmp(x, ROI_temp)), elecs, 'UniformOutput', 0)));


%% Parameter Settings
Name = 'CEWD'; 
curROI = ROI_temp_idx; 
slide_window = 0.020; 
slide_step = 0.010;
dist_measure = 'euclidian';
window_average = 'gaussian';
BL_corr = false;
Cktl_blank_rm = false;


%% Load Data
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


%% Marry Behavioral and EEG Data
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

TrialInfo = TrialInfo(sort_idx,:);
TimeVec1024 = Preproc_Data_1.time;

Data_EEG = [Preproc_Data_2.trial; Preproc_Data_1.trial];

load ArtifactRej
for i = 1:size(trial2rm.(Name),1)
    if(~strcmp(trial2rm.(Name){i,2},'all'))
        Data_EEG(trial2rm.(Name){i,2},trial2rm.(Name){i,1},:) = NaN;
    end
end
curROI(ismember(curROI, cell2mat(trial2rm.(Name)(strcmp(trial2rm.(Name)(:,2),'all'),1)))) = [];

% Baseline Correction
BL_wind = [-0.26 -0.05];
[~,BL_wind_idx] = min(abs(bsxfun(@minus, [TimeVec1024; TimeVec1024], BL_wind')),[],2);
BL_mean = mean(Data_EEG(:,:,BL_wind_idx(1):BL_wind_idx(2)),3);
for ch = 1:128
    na_idx = find(~sum(isnan(BL_mean(:,ch)),2));
    X = BL_mean(na_idx,ch); %[ones(size(BL_mean,1),1) BL_mean(:,ch)];
    Data_EEG(na_idx,ch,:) = squeeze(Data_EEG(na_idx,ch,:)) - X*((X'*X)\X'*squeeze(Data_EEG(na_idx,ch,:)));
end

% Fill in NaN Trials
for i = 1:length(Nan_trial_idx)
    Data_EEG = [Data_EEG(1:(Nan_trial_idx(i)-1),:,:); nan(1,size(Data_EEG,2),size(Data_EEG,3)); Data_EEG(Nan_trial_idx(i):end,:,:)];
end

% Sort Data
Data_EEG = Data_EEG(sort_idx,:,:);

%% Test if align
trl_comp = nan(size(TrialInfo,1),1);
for i = 1:size(TrialInfo,1)
    if(TrialInfo{i,1} > 0)
        if(strcmp(TrialInfo{i,6},'Drawing'))
            trl_comp(i,1) = nansum(nansum(squeeze(Preproc_Data_2.trial(cell2mat(Preproc_Data_2.trialinfo(:,2)) == TrialInfo{i,2},:,:) - Data_EEG(i,:,:)),1),2);
        else
            trl_comp(i,1) = nansum(nansum(squeeze(Preproc_Data_1.trial(cell2mat(Preproc_Data_1.trialinfo(:,2)) == TrialInfo{i,2},:,:) - Data_EEG(i,:,:)),1),2);
        end
    end
end


%% Create and Plot EEG

elec = curROI;
BL_wind = [-0.26 -0.05];
[~,BL_wind_idx] = min(abs(bsxfun(@minus, [TimeVec1024; TimeVec1024], BL_wind')),[],2);
BL_mean = nanmean(Data_EEG(:,:,BL_wind_idx(1):BL_wind_idx(2)),3);
na_idx = find(~sum(isnan(BL_mean),2));

Data_EEG_nanrm = Data_EEG(na_idx,:,:);
BL_mean(sum(isnan(BL_mean),2) > 0,:) = [];

Data_EEG_200bl = bsxfun(@minus, Data_EEG_nanrm, BL_mean);

Data_EEG_linbl = zeros(size(Data_EEG_nanrm));
for ch = 1:128
    X = BL_mean(:,ch); %[ones(size(BL_mean,1),1) BL_mean(:,ch)];
    Data_EEG_linbl(:,ch,:) = squeeze(Data_EEG_nanrm(:,ch,:)) - X*((X'*X)\X'*squeeze(Data_EEG_nanrm(:,ch,:)));
end

Data_EEG_res = zeros(size(Data_EEG_nanrm));
Res_Cov = zeros(size(Data_EEG_nanrm,1),size(Data_EEG_nanrm,1), size(Data_EEG_nanrm,2));
for ch = 1:128
    X = BL_mean(:,ch); %[ones(size(BL_mean,1),1) BL_mean(:,ch)];
    Data_EEG_res(:,ch,:) = squeeze(Data_EEG_nanrm(:,ch,:)) - X*((X'*X)\X'*squeeze(Data_EEG_nanrm(:,ch,:)));
    Res_Cov(:,:,ch) = (1/size(Data_EEG_res,3))*squeeze(Data_EEG_res(:,1,:))*squeeze(Data_EEG_res(:,1,:))';
end

pinv(Res_Cov(:,:,3));

% GA ERP
figure
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_nanrm(:,elec,:),1),2))); 
hold on
% GA ERP Baseline Correction 1
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_200bl(:,elec,:),1),2)))
% GA ERP Baseline Correction 2
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_linbl(:,elec,:),1),2)))

% Drawing vs Picture ERP Baseline Correction (-200 to -50 ms)
figure
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_nanrm(cell2mat(TrialInfo(na_idx,5))==1,elec,:),1),2)),'b')
hold on
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_200bl(cell2mat(TrialInfo(na_idx,5))==1,elec,:),1),2)),'r')
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_linbl(cell2mat(TrialInfo(na_idx,5))==1,elec,:),1),2)),'k')
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_nanrm(cell2mat(TrialInfo(na_idx,5))==2,elec,:),1),2)),'.b')
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_200bl(cell2mat(TrialInfo(na_idx,5))==2,elec,:),1),2)),'.r')
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_linbl(cell2mat(TrialInfo(na_idx,5))==2,elec,:),1),2)),'.k')
hold off

% Animate vs Inanimate ERP Baseline Correction (-200 to -50 ms)
figure
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_nanrm(cell2mat(TrialInfo(na_idx,7))==1,elec,:),1),2)),'b')
hold on
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_200bl(cell2mat(TrialInfo(na_idx,7))==1,elec,:),1),2)),'r')
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_linbl(cell2mat(TrialInfo(na_idx,7))==1,elec,:),1),2)),'k')
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_nanrm(cell2mat(TrialInfo(na_idx,7))==2,elec,:),1),2)),'.b')
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_200bl(cell2mat(TrialInfo(na_idx,7))==2,elec,:),1),2)),'.r')
plot(TimeVec1024, squeeze(nanmean(nanmean(Data_EEG_linbl(cell2mat(TrialInfo(na_idx,7))==2,elec,:),1),2)),'.k')
hold off


%% Create RSA Matrix

% Get Important variables
chan_list = curROI; % transpose(1:length(elecs))
slide_window_s = slide_window; %s
slide_step_s = slide_step; %s
if(strcmp(window_average,'gaussian'))
    average_kern = @(x, dim, gsize) sum(bsxfun(@times, x, (gausswin(gsize)./sum(gausswin(gsize)))'),dim);
    av = 'gaussian'; % Gaussian Average
else
    average_kern = @(x, dim, gsize) mean(x,dim);
    av = 'uniform';  % Uniform Average
end

samplingrate = 1024; %Hz
slide_window_idx = round(samplingrate*(slide_window_s/2));
slide_step_idx = round(samplingrate*slide_step_s);
TimeVec_idx = 1:slide_step_idx:size(Data_EEG,3);
TimeVec = TimeVec1024(TimeVec_idx);

%n = 21;
%figure;plot(1:n, gausswin(n)'./sum(gausswin(n)));hold on; plot(1:n, (1/n)*ones(1,n));ylim([0 1.2/sum(gausswin(n))]);xlim([1 n]); hold off

nonnan_ind = find(cell2mat(TrialInfo(:,1)) ~= 0);
Sigma = zeros(16,length(curROI),length(curROI));
for c = 1:16
    tmp = zeros(size(Data_EEG,3),length(curROI),length(curROI));
    for t = 1:size(Data_EEG,3)
        dat = Data_EEG(intersect(nonnan_ind, (c-1)*8+(1:8)),curROI,t);
        tmp(t,:,:) = cov1para(dat(~sum(isnan(dat),2),:));
    end
    Sigma(c,:,:) = nanmean(tmp,1);
end
tmp = squeeze(mean(Sigma,1));
sigma_inv = tmp^-0.5;
Data_EEG_corr = zeros(size(Data_EEG,1), length(curROI), size(Data_EEG,3));
for t = 1:size(Data_EEG,3)
    Data_EEG_corr(:,:,t) = Data_EEG(:,curROI,t)*sigma_inv;
end

% Full 128 x 128 Matrix
if(strcmp(dist_measure, 'pearson') || strcmp(dist_measure, 'spearman') || strcmp(dist_measure, 'cosine'))
    RSA_Mat_full = repmat(eye(size(Data_EEG,1)),[1 1 length(TimeVec)]);
else
    RSA_Mat_full = repmat(zeros(size(Data_EEG,1)),[1 1 length(TimeVec)]);
end
MDS_Mat_full = zeros(size(Data_EEG,1),2,length(TimeVec));
cur_trial1 = zeros(size(Data_EEG,2),1);
cur_trial2 = zeros(size(Data_EEG,2),1);
fprintf('Progress:  0.0 %%')
txt = sprintf('Progress:  0.0 %%');
for tp = 1:length(TimeVec) 
    time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, size(Data_EEG,3));
    all_patterns = zeros(size(Data_EEG,1), length(curROI));
    for i = 1:size(Data_EEG,1)
        all_patterns(i,:) = average_kern(squeeze(Data_EEG(i,curROI,time_window)), 2, length(time_window))';
    end
    if(BL_corr)
        na_idx = find(~sum(isnan(all_patterns),2));
        X = mean(BL_mean(na_idx,curROI),2); % [ones(length(na_idx),1) mean(BL_mean(na_idx,curROI),2)];
        Res = all_patterns(na_idx,:) - X*((X'*X)\X'*all_patterns(na_idx,:));
        Cov_Res = (1/length(na_idx))*(Res'*Res);
    end
    if(Cktl_blank_rm)
        mean_pattern = nanmean(all_patterns,1)';
    end
    
    for i = 1:size(Data_EEG,1)-1
        for j = (i+1):size(Data_EEG,1)
            cur_trial1 = average_kern(squeeze(Data_EEG(i,curROI,time_window)), 2, length(time_window));
            cur_trial2 = average_kern(squeeze(Data_EEG(j,curROI,time_window)), 2, length(time_window));
            if(Cktl_blank_rm && (strcmp(dist_measure, 'euclidian') || strcmp(dist_measure, 'mahalanobis')))
                cur_trial1 = cur_trial1 - mean_pattern;
                cur_trial2 = cur_trial2 - mean_pattern;
            end
            if(BL_corr)
                RSA_Mat_full(i,j,tp) = (cur_trial1 - cur_trial2)'/Cov_Res*(cur_trial1 - cur_trial2);
                RSA_Mat_full(j,i,tp) = (cur_trial1 - cur_trial2)'/Cov_Res*(cur_trial1 - cur_trial2);
            else
                RSA_Mat_full(i,j,tp) = dist_meas(cur_trial1, cur_trial2, dist_measure);
                RSA_Mat_full(j,i,tp) = dist_meas(cur_trial1, cur_trial2, dist_measure);
            end
        end
    end
    
    RSA_Mat_full(RSA_Mat_full < 0) = NaN;
    RSA_Mat_full(:,:,tp) = sqrt(RSA_Mat_full(:,:,tp));
    
    if(tp == 1)
        MDS_Mat_full(:,:,tp) = mdscale(RSA_Mat_full(:,:,tp),2,'Start','random');
    else
        MDS_Mat_full(:,:,tp) = mdscale(RSA_Mat_full(:,:,tp),2,'Start',MDS_Mat_full(:,:,tp-1));
        %[U,S,V] = svd(MDS_Mat_av(:,:,tp-1)'*MDS_Mat_av(:,:,tp));
        %r = sign(det(V*U'));
        %R = V*[1 0; 0 r]*U';
        %MDS_Mat_av(:,:,tp) = MDS_Mat_av(:,:,tp)*R;
    end
    
    fprintf(repmat('\b',1,length(txt)))
    fprintf('Progress: %3.2f %%',(tp/length(TimeVec))*100)
    txt = sprintf('Progress: %3.2f %%',(tp/length(TimeVec))*100);
end
fprintf('\n')

figure
time_plot = -0.2:0.1:1.5;
for sbp = 1:18
    [~,t_idx] = min(abs(TimeVec - time_plot(sbp)));
    subplot(3,6,sbp)
    imagesc(squeeze(RSA_Mat_full(:,:,t_idx)))
    colorbar; title(sprintf('%4.0f ms', TimeVec(t_idx)*1000))
    axis square
end

figure
[~,t_idx] = min(abs(TimeVec - 0.3));
imagesc(squeeze(RSA_Mat_full(:,:,t_idx))); colorbar
title(sprintf('RSA Matrix at %4.0f ms', TimeVec(t_idx)*1000))
axis square

figure
time_plot = -0.2:0.1:1.5;
for sbp = 1:18
    [~,t_idx] = min(abs(TimeVec - time_plot(sbp)));
    subplot(3,6,sbp)
    gscatter(MDS_Mat_full(:,1,t_idx), MDS_Mat_full(:,2,t_idx), cell2mat(TrialInfo(:,7)))
    title(sprintf('%4.0f ms', TimeVec(t_idx)*1000))
    ylim([-150 150]); xlim([-150 150])
    axis square
end


% Average over Examples: 16 x 16 Matrix
if(strcmp(dist_measure, 'pearson') || strcmp(dist_measure, 'spearman') || strcmp(dist_measure, 'cosine'))
    RSA_Mat_av = repmat(eye(16),[1 1 length(TimeVec)]);
else
    RSA_Mat_av = repmat(zeros(16),[1 1 length(TimeVec)]);
end
MDS_Mat_av = zeros(16,3,length(TimeVec));
%Data = Data_EEG(:,curROI,:);
Data = Data_EEG_corr;
ind1 = [];
ind2 = [];
cur_trial1 = zeros(8,length(chan_list));
cur_trial2 = zeros(8,length(chan_list));
cur_BL1 = zeros(8,length(chan_list));
cur_BL2 = zeros(8,length(chan_list));
BL_mean = nanmean(Data_EEG(:,curROI,BL_wind_idx(1):BL_wind_idx(2)),3);
fprintf('Progress:  0.0 %%')
txt = sprintf('Progress:  0.0 %%');
for tp = 1:length(TimeVec) 
    time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, size(Data,3));
    all_patterns = zeros(size(Data,1), length(curROI));
    for i = 1:size(Data,1)
        all_patterns(i,:) = average_kern(squeeze(Data(i,:,time_window)), 2, length(time_window))';
    end
    if(Cktl_blank_rm)
        mean_pattern = nanmean(all_patterns,1)';
    end
    
    dat1 = zeros(8,length(chan_list),length(time_window));
    dat2 = zeros(8,length(chan_list),length(time_window));
    ct_x = 1; ct_y = 1;
    for i = 1:8:size(Data,1)-1
        for j = i:8:size(Data,1)
            if(i ~= j)
                for h = 0:7
                    cur_trial1(h+1,:) = average_kern(squeeze(Data(i+h,:,time_window)), 2, length(time_window))';
                    cur_trial2(h+1,:) = average_kern(squeeze(Data(j+h,:,time_window)), 2, length(time_window))';
                end
                cur_trial1(sum(isnan(cur_trial1),2) > 0,:) = []; cur_trial2(sum(isnan(cur_trial2),2) > 0,:) = []; Cur_trial = [cur_trial1; cur_trial2];

    %             % Residual Variance
    %             Res = zeros(size(Dat,2));
    %             for n = 1:size(Dat,1)
    %                 X = [ones(length(curROI),1) cur_BL(n,:)'];
    %                 y = squeeze(Dat(n,:,:));
    %                 R = y - X*((X'/squeeze(Sigma(Ind(n),:,:))*X)\(X'/squeeze(Sigma(Ind(n),:,:))*y));
    %                 Res = Res + R*R';
    %             end
    %             Res = Res./(size(Dat,1)*(length(time_window) - size(X,2)));

                % Create all differences
                poss_folds = [kron(1:size(cur_trial1,1),ones(1,size(cur_trial2,1)))' kron(ones(1,size(cur_trial1,1)),1:size(cur_trial2,1))'];
                diff_crossv = zeros(size(cur_trial1,2), size(poss_folds,1));
                for df = 1:size(poss_folds,1)
                    diff_crossv(:,df) = cur_trial1(poss_folds(df,1),:)' - cur_trial2(poss_folds(df,2),:)';
                end

                % Crossvalidation Folds
                dist_crossv = zeros(size(diff_crossv,1),2);
                for fold = 1:size(diff_crossv,2)
                    dist_crossv(fold,1) = diff_crossv(:,fold)'*mean(diff_crossv(:,[1:fold-1 fold+1:end]),2);
                    %dist_crossv(fold,2) = diff_crossv(:,fold)'*sigma_inv*mean(diff_crossv(:,[1:fold-1 fold+1:end]),2);
                end
                
                % SVM and LDA
                cfg = [];
                cfg.metric      = 'acc';
                cfg.cv          = 'leaveout';
                cfg.repeat      = 2;
                cfg.classifier  = 'lda';
                cfg.feedback    = 0;
                lda_acc_class   = mv_classify(cfg, Cur_trial, [ones(1,size(cur_trial1,1)) 2*ones(1,size(cur_trial2,1))]');
                
                cfg.classifier  = 'svm';
                svm_acc_class   = mv_classify(cfg, Cur_trial, [ones(1,size(cur_trial1,1)) 2*ones(1,size(cur_trial2,1))]');
                
                RSA_Mat_av(ct_x,ct_y,tp) = lda_acc_class; %mean(dist_crossv(:,3));
                RSA_Mat_av(ct_y,ct_x,tp) = lda_acc_class; %mean(dist_crossv(:,3)); 
            end
            ct_y = ct_y + 1;
        end
        ct_x = ct_x + 1;
        ct_y = ct_x;
    end
    RSA_Mat_av(RSA_Mat_av < 0) = NaN;
    %RSA_Mat_av(:,:,tp) = sqrt(RSA_Mat_av(:,:,tp));
    
    if(tp == 1)
        MDS_Mat_av(:,:,tp) = mdscale(RSA_Mat_av(:,:,tp),3,'Start','random');
    else
        MDS_Mat_av(:,:,tp) = mdscale(RSA_Mat_av(:,:,tp),3,'Start',MDS_Mat_av(:,:,tp-1));
        %[U,S,V] = svd(MDS_Mat_av(:,:,tp-1)'*MDS_Mat_av(:,:,tp));
        %r = sign(det(V*U'));
        %R = V*[1 0; 0 r]*U';
        %MDS_Mat_av(:,:,tp) = MDS_Mat_av(:,:,tp)*R;
    end
    
    
    fprintf(repmat('\b',1,length(txt)))
    fprintf('Progress: %3.2f %%',(tp/length(TimeVec))*100)
    txt = sprintf('Progress: %3.2f %%',(tp/length(TimeVec))*100);
end
fprintf('\n')

figure
time_plot = -0.2:0.1:1.5;
for sbp = 1:18
    [~,t_idx] = min(abs(TimeVec - time_plot(sbp)));
    subplot(3,6,sbp)
    imagesc(squeeze(RSA_Mat_av(:,:,t_idx))); ax = gca; %set(ax,'clim',[0 20])
    colorbar; title(sprintf('%4.0f ms', TimeVec(t_idx)*1000))
    axis square
end

figure
[~,t_idx] = min(abs(TimeVec - 0.3));
imagesc(squeeze(RSA_Mat_av(:,:,t_idx))); colorbar
title(sprintf('RSA Matrix at %4.0f ms', TimeVec(t_idx)*1000))
axis square

figure
time_plot = -0.2:0.1:1.5;
for sbp = 1:18
    [~,t_idx] = min(abs(TimeVec - time_plot(sbp)));
    subplot(3,6,sbp)
    gscatter(MDS_Mat_av(:,1,t_idx), MDS_Mat_av(:,2,t_idx), [ones(4,1); 2*ones(4,1); ones(4,1); 2*ones(4,1)])   % [ones(8,1); 2*ones(8,1)]
    title(sprintf('%4.0f ms', TimeVec(t_idx)*1000))
    %ylim([-1000 1000]); xlim([-1000 1000])
    axis square
end

figure
plot3(squeeze(mean(MDS_Mat_av([1:4 9:12],1,:),1)), squeeze(mean(MDS_Mat_av([1:4 9:12],2,:),1)), TimeVec)
hold on
plot3(squeeze(mean(MDS_Mat_av([5:8 13:16],1,:),1)), squeeze(mean(MDS_Mat_av([5:8 13:16],2,:),1)), TimeVec)




%% t-SNE

tSNE_Mat = repmat(eye(size(Data_EEG,1)),[1 1 length(TimeVecNew)]);
cur_trial = zeros(size(Data_EEG,1),length(chan_list));
fprintf('Progress:  0.0 %%')
txt = sprintf('Progress:  0.0 %%');
for tp = 1:length(TimeVecNew) 
    time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, size(Data_EEG,3));
    for i = 1:size(Data_EEG,1)
        cur_trial(i,:) = average_kern(squeeze(Data_EEG(i,chan_list,time_window)), 2, length(time_window))';
    end
    tSNE_Mat = tsne(cur_trial);
    fprintf(repmat('\b',1,length(txt)))
    fprintf('Progress: %3.2f %%',(tp/length(TimeVecNew))*100)
    txt = sprintf('Progress: %3.2f %%',(tp/length(TimeVecNew))*100);
end
fprintf('\n')


%% Hypotheses Matrix

Perceptual_Mat_full = zeros(size(TrialInfo,1));
Semantic_Mat_full = zeros(size(TrialInfo,1));
for i = 1:size(TrialInfo,1)-1
    for j = (i+1):size(TrialInfo,1)
        if(j ~= size(TrialInfo,1) - (i - 1))
            if(TrialInfo{i,5} == 1 && TrialInfo{j,5} == 1)
                Perceptual_Mat_full(i,j) = 1;
            elseif(TrialInfo{i,5} == 2 && TrialInfo{j,5} == 2)
                Perceptual_Mat_full(i,j) = 2;
            else
                Perceptual_Mat_full(i,j) = -1;
            end

            if(TrialInfo{i,7} == 1 && TrialInfo{j,7} == 1)
                Semantic_Mat_full(i,j) = 1;
            elseif(TrialInfo{i,7} == 2 && TrialInfo{j,7} == 2)
                Semantic_Mat_full(i,j) = 2;
            else
                Semantic_Mat_full(i,j) = -1;
            end
        end
    end
end

figure
subplot(1,2,1)
imagesc(Perceptual_Mat_full); title('Perceptual Hypothesis Matrix')
axis square
subplot(1,2,2)
imagesc(Semantic_Mat_full); title('Semantic Hypothesis Matrix')
axis square


Perceptual_Mat_av = zeros(16);
Semantic_Mat_av = zeros(16);
for i = 1:16-1
    for j = (i+1):16
        if(j ~= 16 - (i - 1))
            if(TrialInfo{8*i,5} == 1 && TrialInfo{8*j,5} == 1)
                Perceptual_Mat_av(i,j) = 1;
            elseif(TrialInfo{8*i,5} == 2 && TrialInfo{8*j,5} == 2)
                Perceptual_Mat_av(i,j) = 2;
            else
                Perceptual_Mat_av(i,j) = -1;
            end

            if(TrialInfo{8*i,7} == 1 && TrialInfo{8*j,7} == 1)
                Semantic_Mat_av(i,j) = 1;
            elseif(TrialInfo{8*i,7} == 2 && TrialInfo{8*j,7} == 2)
                Semantic_Mat_av(i,j) = 2;
            else
                Semantic_Mat_av(i,j) = -1;
            end
        end
    end
end

figure
subplot(1,2,1)
imagesc(Perceptual_Mat_av); title('Perceptual Hypothesis Matrix')
axis square
subplot(1,2,2)
imagesc(Semantic_Mat_av); title('Semantic Hypothesis Matrix')
axis square



%% Perceptual Dimension

RSA_Mat = RSA_Mat_av;
Perceptual_Mat = Perceptual_Mat_av;

cur_data = zeros(size(RSA_Mat,1));
RSA_Perceptual_drawing = zeros(1,length(TimeVec));
RSA_Perceptual_picture = zeros(1,length(TimeVec));
RSA_Perceptual_within  = zeros(1,length(TimeVec));
RSA_Perceptual_between = zeros(1,length(TimeVec));
for tp = 1:length(TimeVec)
    cur_data = RSA_Mat(:,:,tp);
    RSA_Perceptual_drawing(1,tp) = nanmean(cur_data(Perceptual_Mat == 1));
    RSA_Perceptual_picture(1,tp) = nanmean(cur_data(Perceptual_Mat == 2));
    RSA_Perceptual_within(1,tp) = nanmean(cur_data(Perceptual_Mat > 0));
    RSA_Perceptual_between(1,tp) = nanmean(cur_data(Perceptual_Mat < 0));
end

figure
plot(TimeVec, RSA_Perceptual_drawing)
hold on
plot(TimeVec, RSA_Perceptual_picture)
hold off
ylabel('Corr (z)'); xlabel('Time (s)'); title('Perceptual Dim Drawing vs. Picture')
xlim([-0.1 1.5]); legend('drawing','picture')

figure
plot(TimeVec, RSA_Perceptual_within)
hold on
plot(TimeVec, RSA_Perceptual_between)
hold off
ylabel('Corr (z)'); xlabel('Time (s)'); title('Perceptual Dim Within vs. Between Cat')
xlim([-0.1 1.5]); legend('perceptual within','perceptual between')


%% Semantic Dimension

RSA_Mat = RSA_Mat_av;
Semantic_Mat = Semantic_Mat_av;

cur_data = zeros(size(RSA_Mat,1));
RSA_Semantic_animate = zeros(1,length(TimeVec));
RSA_Semantic_inanimate = zeros(1,length(TimeVec));
RSA_Semantic_within  = zeros(1,length(TimeVec));
RSA_Semantic_between = zeros(1,length(TimeVec));
for tp = 1:length(TimeVec)
    cur_data = RSA_Mat(:,:,tp);
    RSA_Semantic_animate(1,tp) = mean(cur_data(Semantic_Mat == 1));
    RSA_Semantic_inanimate(1,tp) = mean(cur_data(Semantic_Mat == 2));
    RSA_Semantic_within(1,tp) = mean(cur_data(Semantic_Mat > 0));
    RSA_Semantic_between(1,tp) = mean(cur_data(Semantic_Mat < 0));
end

figure
plot(TimeVec, RSA_Semantic_animate)
hold on
plot(TimeVec, RSA_Semantic_inanimate)
hold off
ylabel('Corr (z)'); xlabel('Time (s)'); title('Semantic Dim Animate vs. Inanimate')
xlim([-0.1 1.5]); legend('Animate','Inanimate')

figure
plot(TimeVec, RSA_Semantic_within)
hold on
plot(TimeVec, RSA_Semantic_between)
hold off
ylabel('Corr (z)'); xlabel('Time (s)'); title('Semantic Dim Within vs. Between Cat')
xlim([-0.1 1.5]); legend('semantic within','semantic between')


%% Perceptual vs Semantic

figure
plot(TimeVec, RSA_Perceptual_within)
hold on
plot(TimeVec, RSA_Semantic_within)
hold off
ylabel('Corr (z)'); xlabel('Time (s)'); title('Perceptual vs. Semantic')
xlim([-0.1 1.5]); legend('Perceptual','Semantic')


