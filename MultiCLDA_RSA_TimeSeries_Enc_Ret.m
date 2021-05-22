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



% Plot RSA Difference Time series

r1 = 1; r2 = 2;
c1 = 1; c2 = 2;

Data1 = MultiCLDA_RSA_Enc_Data.(ROI{r1}).Data;
Data2 = MultiCLDA_RSA_Enc_Data.(ROI{r2}).Data;
TimeVec = MultiCLDA_RSA_Enc_Data.TimeEnc;
Cond_idx = MultiCLDA_RSA_Enc_Data.Cond_idx;

% cfg = [];
% cfg.nPerms = 1000;
% cfg.thresh_pval = 0.05;
% cfg.mcc_cluster_pval = 0.05;
% cfg.TimeVec = RSA_Data.TimeVec;
% cfg.Hyp_Mat = Hyp_Mat_per;
% cfg.matshuff = false;
% cfg.twoside = true;
% Results1 = rsa_perm(cfg, RSA_Data.(ROI{r1}).red16_Data);
% 
% cfg.Hyp_Mat = Hyp_Mat_sem;
% Results2 = rsa_perm(cfg, RSA_Data.(ROI{r2}).red16_Data);

figure('Pos', [325 510 650 402]);
plot([]); hold on
dat1 = mean(Data1(:,:,Cond_idx.([Cat{c},'_WI_BT_idx'])(:,1)),3) - mean(Data1(:,:,Cond_idx.([Cat{c},'_WI_BT_idx'])(:,2)),3);
dat2 = mean(Data2(:,:,Cond_idx.([Cat{c},'_WI_BT_idx'])(:,1)),3) -mean(Data2(:,:,Cond_idx.([Cat{c},'_WI_BT_idx'])(:,2)),3);
SEM1 = nanstd(dat1,0,1)./sqrt(size(dat1,1)); SEM2 = nanstd(dat2,0,1)./sqrt(size(dat2,1));
% SEM1 = zeros(1,length(TimeVec)); SEM2 = zeros(1,length(TimeVec));
% for tp = 1:length(TimeVec)
%     [~,~,~,stats] = ttest(RSA_Time.(ROI{r1}).(Cat{c1}).(Dat_names{dt(2)})(:,tp), RSA_Time.(ROI{r1}).(Cat{c1}).(Dat_names{dt(1)})(:,tp));
%     SEM1(1,tp) = tinv(.95,stats.df)*(stats.sd/sqrt(size(dat1,1)));
%     [~,~,~,stats] = ttest(RSA_Time.(ROI{r2}).(Cat{c2}).(Dat_names{dt(2)})(:,tp), RSA_Time.(ROI{r2}).(Cat{c2}).(Dat_names{dt(1)})(:,tp));
%     SEM2(1,tp) = tinv(.95,stats.df)*(stats.sd/sqrt(size(dat2,1)));
% end
fill([TimeVec fliplr(TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
xlim([-0.2 1.5]); %ylim([0.35 0.6]);
plot([-0.2 1.5],[0 0],'--k','linewidth',1)
h1 = plot(TimeVec, nanmean(dat1,1),'b','linewidth',2);
h2 = plot(TimeVec, nanmean(dat2,1),'r','linewidth',2);
ylabel(['Diff ',RSA_Data.meas16{msr}]); xlabel('Time (s)'); %title([Cat{c},': ',Dat_names{dt(1)},' vs ',Dat_names{dt(2)}])
lg = legend([h1 h2], {[ROI{r1},' ',Cat{c1},' ',Dat_names1{dt1(2)},' - ',Dat_names1{dt1(1)}],[ROI{r2},' ',Cat{c2},' ',Dat_names2{dt2(2)},' - ',Dat_names2{dt2(1)}]}); legend boxoff; set(lg,'FontSize',14)
box off;
set(gca,'linewidth',2.5,'FontSize',14,'xlim',[-0.2 1])
sign_mcc_clust_1 = Results1.zmapthresh;
sign_mcc_clust_1(sign_mcc_clust_1 > 0) = min(get(gca,'ylim'))*0.8;
plot(TimeVec,sign_mcc_clust_1,'bo','MarkerFaceColor','b')
sign_mcc_clust_2 = Results2.zmapthresh;
sign_mcc_clust_2(sign_mcc_clust_2 > 0) = min(get(gca,'ylim'))*0.9;
plot(TimeVec,sign_mcc_clust_2,'ro','MarkerFaceColor','r')
hold off



%% Plot Sensor Weights

r = 3;
elec_idx = ROI_idx{r};

Sensor_Data = zeros(size(MultiCLDA_RSA_Enc_Data.(ROI{r}).Sensor_weights,1), length(ROI_all_idx), size(MultiCLDA_RSA_Enc_Data.(ROI{r}).Sensor_weights,2));
Sensor_Data(:,elec_idx,:) = permute(MultiCLDA_RSA_Enc_Data.(ROI{r}).Sensor_weights,[1 3 2]);

GA_Sensor_Data = [];
GA_Sensor_Data.avg = squeeze(mean(Sensor_Data,1)); 
GA_Sensor_Data.var = squeeze(std(Sensor_Data,0,1));  
GA_Sensor_Data.dof = size(Sensor_Data,1).*ones(size(Sensor_Data,2),size(Sensor_Data,3)); 
GA_Sensor_Data.time = MultiCLDA_RSA_Enc_Data.TimeEnc; 
GA_Sensor_Data.label = elecs; 
GA_Sensor_Data.dimord = 'chan_time';

GA_Sensor_Data.mask = zeros(size(GA_Sensor_Data.avg));
GA_Sensor_Data.mask(elec_idx,:) = 1;

cfg = [];
cfg.xlim             = [-0.2:0.1:1];
cfg.maskparameter    = 'mask'; 
%cfg.colorbar         = 'yes';
cfg.layout           = 'biosemi128.lay';
cfg.highlight        = 'on';
cfg.highlightchannel = elec_idx;
cfg.highlightsymbol  = 's';
cfg.highlightcolor   = [0 0 0];
cfg.comment          = 'no';

ft_topoplotER(cfg, GA_Sensor_Data);







