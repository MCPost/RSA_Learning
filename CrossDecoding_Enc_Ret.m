%% CrossDecoding Encoding and Retrieval

% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs

% Hypotheses Matrices
Hypothesis_Matrix

%% Load Data 
%% Subject Names
load('RSA_Data_Enc','Subj_names')

%% Measure
msr = 1;

%% Matrix Indices to choose
UTmat = triu(ones(16)).*~eye(16);
UTM_idx = find(UTmat(:) == 1);
SensUTmat = UTmat(1:end-1,:);
SensUTM_idx = find(SensUTmat(:) == 1);


%% Cross-correlate RSA time courses for each subject

% cfg Struct
cfg = [];
cfg.subs             = [];
cfg.slide_window     = 0.080;
cfg.slide_step       = 0.010;
cfg.fwhm             = 0.040;
cfg.window_average   = 'gaussian';
cfg.MNN              = true;
cfg.Hyp_perceptual   = {Perceptual_Mat_full Perceptual_Mat_red16}; %{triu(ones(128)).*~eye(128) triu(ones(16)).*~eye(16)}; %{Perceptual_Mat_full Perceptual_Mat_red16};
cfg.Hyp_semantic     = {Semantic_Mat_full   Semantic_Mat_red16};
cfg.ROI              = {};
cfg.pair             = '';
cfg.msr              = [];
cfg.timewind1        = [-0.2 1.5];
cfg.timewind2        = [-2.5 0.2];
cfg.only16           = true;
cfg.permtest         = false;
cfg.n_perms          = 1000;
cfg.thresh_pval      = 0.05;
cfg.mcc_cluster_pval = 0.05;
cfg.ts_os_fac        = 1;
cfg.matshuffle       = true;

jobHandles = {};


% Encoding to Retrieval, Between ROI, No PermTest
cfg.permtest         = true; %false;
cfg.timewind2        = [-2.2 -0.8];
cfg.ROI              = {'OCC','TMP';...
                        %'OCC','FRT';...
                        'OCC','PRT';...
                        'TMP','OCC';...
                        %'TMP','FRT';...
                        'TMP','PRT';...
                        %'FRT','OCC';...
                        %'FRT','TMP';...
                        %'FRT','PRT';...
                        'PRT','OCC';...
                        'PRT','TMP';...
                        'PRT','FRT'};
pair = 'Enc - Ret';
save_name = 'Crossdecoding_Matfiles/CrossDec_RSA_Enc_Ret_btROI_perm_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1
    jobHandles{end_ind+msr} = batch(c, @CrossDecoding_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
end

jobHandles{1}
diary(jobHandles{1})
%delete(jobHandles{1})


% Encoding to Encoding, Within ROI, No PermTest
cfg.permtest         = false;
cfg.timewind2        = [-0.2 1.5];
cfg.ROI              = {'OCC';...
                        'TMP';...
                        'FRT';...
                        'PRT'};
pair = 'Enc - Enc';
save_name = 'Crossdecoding_Matfiles/CrossDec_RSA_Enc_Enc_wiROI_noperm_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1
    jobHandles{end_ind+msr} = batch(c, @CrossDecoding_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
end


% Encoding to Encoding, Between ROI, No PermTest
cfg.permtest         = false;
cfg.timewind2        = [-0.2 1.5];
cfg.ROI              = {'OCC','TMP';'OCC','FRT';'OCC','PRT';...
                        'TMP','FRT';'TMP','PRT';...
                        'FRT','PRT'};
pair = 'Enc - Enc';
save_name = 'Crossdecoding_Matfiles/CrossDec_RSA_Enc_Enc_btROI_noperm_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1
    jobHandles{end_ind+msr} = batch(c, @CrossDecoding_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
end

jobHandles{1}
diary(jobHandles{1})
%delete(jobHandles{1})


% Retrieval to Retrieval, Within ROI, No PermTest
cfg.permtest         = false;
cfg.timewind2        = [-2.5 0.2];
cfg.ROI              = {'OCC';...
                        'TMP';...
                        'FRT';...
                        'PRT'};
pair = 'Ret - Ret';
save_name = 'Crossdecoding_Matfiles/CrossDec_RSA_Ret_Ret_wiROI_noperm_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1
    jobHandles{end_ind+msr} = batch(c, @CrossDecoding_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
end


% Retrieval to Retrieval, Between ROI, No PermTest
cfg.permtest         = false;
cfg.timewind2        = [-2.5 0.2];
cfg.ROI              = {'OCC','TMP';'OCC','FRT';'OCC','PRT';...
                        'TMP','FRT';'TMP','PRT';...
                        'FRT','PRT'};
pair = 'Ret - Ret';
save_name = 'Crossdecoding_Matfiles/CrossDec_RSA_Ret_Ret_btROI_noperm_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1
    jobHandles{end_ind+msr} = batch(c, @CrossDecoding_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
end

jobHandles{1}
diary(jobHandles{1})
%delete(jobHandles{1})


% Encoding to Retrieval, Within ROI, PermTest
cfg.permtest         = false; %true;
cfg.timewind2        = [-2.5 0.2];
%cfg.ROI              = {'Wholehead'};
cfg.ROI              = {'OCC';...
                        'TMP';...
                        'FRT';...
                        'PRT'};

cfg.pair = 'Enc - Ret';
save_name = 'Crossdecoding_Matfiles/CrossDec_Enc_Ret_wiROI_noperm_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1
    jobHandles{end_ind+msr} = batch(c, @CrossDecoding_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
end

jobHandles{4}
diary(jobHandles{4})





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

time_wind_enc = [-0.1  1.2];
time_wind_ret = [-1.8  -0.6];

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






