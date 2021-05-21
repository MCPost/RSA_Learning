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

% Load Enc Data and Ret Data
load('MultiCLDA_RSA_Enc')
load('MultiCLDA_RSA_Ret')


% Matrix Indices to choose
UTmat = triu(ones(16)).*~eye(16);
UTM_idx = find(UTmat(:) == 1);

% Reshape to RDM
hlpmat = zeros(16);
RDM_Enc = []; RDM_Ret = [];
for r = 1:length(ROI)
    RDM_Enc.(ROI{r}) = zeros(size(MultiCLDA_RSA_Enc_Data.(ROI{r}).Data,1),size(MultiCLDA_RSA_Enc_Data.(ROI{r}).Data,1),16,16);
    for s = 1:size(MultiCLDA_RSA_Enc_Data.(ROI{r}).Data,1)
        for t = 1:size(MultiCLDA_RSA_Enc_Data.(ROI{r}).Data,1)
            hlpmat(UTM_idx) = MultiCLDA_RSA_Enc_Data.(ROI{r}).Data(s,t,:);
            RDM_Enc.(ROI{r})(s,t,:,:) = hlpmat;
        end
    end
    RDM_Ret.(ROI{r}) = zeros(size(MultiCLDA_RSA_Ret_Data.(ROI{r}).Data,1),size(MultiCLDA_RSA_Ret_Data.(ROI{r}).Data,1),16,16);
    for s = 1:size(MultiCLDA_RSA_Ret_Data.(ROI{r}).Data,1)
        for t = 1:size(MultiCLDA_RSA_Ret_Data.(ROI{r}).Data,1)
            hlpmat(UTM_idx) = MultiCLDA_RSA_Ret_Data.(ROI{r}).Data(s,t,:);
            RDM_Ret.(ROI{r})(s,t,:,:) = hlpmat;
        end
    end
end


% Plot RSA Time series

r = 1;
Cat = {'Perc','Sem'}; c = 1;

Data = MultiCLDA_RSA_Enc_Data.(ROI{r}).Data;
TimeVec = MultiCLDA_RSA_Enc_Data.TimeEnc;
Cond_idx = MultiCLDA_RSA_Enc_Data.Cond_idx;

figure('Pos', [325 510 650 402]);
plot([]); hold on
dat1 = mean(mean(Data(:,:,Cond_idx.([Cat{c},'_WI_BT_idx'])(:,1)),3),1);
dat2 = mean(mean(Data(:,:,Cond_idx.([Cat{c},'_WI_BT_idx'])(:,2)),3),1);
SEM1 = nanstd(mean(Data(:,:,Cond_idx.([Cat{c},'_WI_BT_idx'])(:,1)),3),0,1)./sqrt(size(Data,1));
SEM2 = nanstd(mean(Data(:,:,Cond_idx.([Cat{c},'_WI_BT_idx'])(:,2)),3),0,1)./sqrt(size(Data,1));
fill([TimeVec fliplr(TimeVec)],[dat1 fliplr(dat1 + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[dat1 fliplr(dat1 - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[dat2 fliplr(dat2 + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[dat2 fliplr(dat2 - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
h1 = plot(TimeVec, dat1,'b','linewidth',2);
h2 = plot(TimeVec, dat2,'r','linewidth',2);
hold off
ylabel('LDA'); xlabel('Time (s)')
xlim([-0.2 1.5]); %lg = legend([h1 h2], {Dat_names{dt(1)},Dat_names{dt(2)}}); legend boxoff; set(lg,'FontSize',14)
box off; %ylim([0.42 0.6])
set(gca,'linewidth',2.5,'FontSize',14)



