%% CrossDecoding Encoding and Retrieval


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs 
ROI = {'OCC';'TMP';'FRT';'PRT'};
ROI_idx = {ROI_occ_idx; ROI_temp_idx; ROI_front_idx; ROI_pari_idx};

% Subject Names
load('RSA_Data_Enc','Subj_names')

% Load Data
tmp_struct = load('CrossDec_Enc_Ret');
tmp_fnames = fieldnames(tmp_struct);

time_wind_enc = [-0.1  1.5];
time_wind_ret = [-2.5  0.1];

Cat = {'Perc','Sem'};
ROI = tmp_struct.(tmp_fnames{1}).ROI;
r1 = 1; r2 = 2;
c1 = 1; c2 = 2;

TimeX = tmp_struct.(tmp_fnames{1}).TimeX;
time_idx_enc = dsearchn(TimeX',time_wind_enc')';
TimeY = tmp_struct.(tmp_fnames{1}).TimeY;
time_idx_ret = dsearchn(TimeY',time_wind_ret')';
TimeX = TimeX(time_idx_enc(1):time_idx_enc(end));
TimeY = TimeY(time_idx_ret(1):time_idx_ret(end));
CurDat1 = tmp_struct.(tmp_fnames{1}).(ROI{r1}).(['Data_',Cat{c1}])(:,time_idx_ret(1):time_idx_ret(end),time_idx_enc(1):time_idx_enc(end));
CurDat2 = tmp_struct.(tmp_fnames{1}).(ROI{r2}).(['Data_',Cat{c2}])(:,time_idx_ret(1):time_idx_ret(end),time_idx_enc(1):time_idx_enc(end));



%% Compare Perceptual to Semantic Dimension

enc_lim = time_wind_enc;
ret_lim = time_wind_ret;


cfg = [];
cfg.nPerms = 1000;
cfg.thresh_pval = 0.05;
cfg.mcc_cluster_pval = 0.05;
cfg.twoside = false;
Results5 = PermTestRNDCP2D(cfg, CurDat1);
Results6 = PermTestRNDCP2D(cfg, CurDat2);


figure('Pos',[189 137 1531 782])
subplot(1,2,1)
curdat = squeeze(mean(CurDat1,1));
contourf(TimeX, TimeY, curdat, 40,'linestyle','none'); colorbar
caxis(prctile(curdat(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim); 
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
if(Results5.H == 1)
    zmapthresh = Results5.zmapthresh; zmapthresh(isnan(zmapthresh)) = 0; zmapthresh(zmapthresh < 0) = 0; zmapthresh(zmapthresh ~= 0) = 1;
    contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
end
hold off
subplot(1,2,2)
curdat = squeeze(mean(CurDat2,1));
contourf(TimeX, TimeY, curdat, 40,'linestyle','none'); colorbar
caxis(prctile(curdat(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim); 
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
if(Results6.H == 1)
    zmapthresh = Results6.zmapthresh; zmapthresh(isnan(zmapthresh)) = 0; zmapthresh(zmapthresh < 0) = 0; zmapthresh(zmapthresh ~= 0) = 1;
    contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
end
hold off






