%% Cross-comparison of RSA time courses

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
cfg.subs             = Subj_names;
cfg.slide_window     = 0.080;
cfg.slide_step       = 0.010;
cfg.fwhm             = 0.040;
cfg.window_average   = 'gaussian';
cfg.Hyp_perceptual   = {Perceptual_Mat_full Perceptual_Mat_red16}; %{triu(ones(128)).*~eye(128) triu(ones(16)).*~eye(16)}; %{Perceptual_Mat_full Perceptual_Mat_red16};
cfg.Hyp_semantic     = {Semantic_Mat_full   Semantic_Mat_red16};
cfg.ROI              = {};
cfg.timewind1        = [-0.2 1.5];
cfg.timewind2        = [-2.5 0.2];
cfg.only16           = true;
cfg.permtest         = false;
cfg.n_perms          = 1000;
cfg.thresh_pval      = 0.05;
cfg.mcc_cluster_pval = 0.05;
cfg.ts_os_fac        = 1;
cfg.matshuffle       = true;
cfg.studentized      = false;

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
save_name = 'Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Ret_btROI_perm_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1:3
    jobHandles{end_ind+msr} = batch(c, @Crosscomparison_RSA_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
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
save_name = 'Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Enc_wiROI_noperm_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1:3
    jobHandles{end_ind+msr} = batch(c, @Crosscomparison_RSA_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
end


% Encoding to Encoding, Between ROI, No PermTest
cfg.permtest         = false;
cfg.timewind2        = [-0.2 1.5];
cfg.ROI              = {'OCC','TMP';'OCC','FRT';'OCC','PRT';...
                        'TMP','FRT';'TMP','PRT';...
                        'FRT','PRT'};
pair = 'Enc - Enc';
save_name = 'Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Enc_btROI_noperm_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1:3
    jobHandles{end_ind+msr} = batch(c, @Crosscomparison_RSA_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
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
save_name = 'Crosscompare_RSA_Matfiles/CrossComp_RSA_Ret_Ret_wiROI_noperm_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1:3
    jobHandles{end_ind+msr} = batch(c, @Crosscomparison_RSA_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
end


% Retrieval to Retrieval, Between ROI, No PermTest
cfg.permtest         = false;
cfg.timewind2        = [-2.5 0.2];
cfg.ROI              = {'OCC','TMP';'OCC','FRT';'OCC','PRT';...
                        'TMP','FRT';'TMP','PRT';...
                        'FRT','PRT'};
pair = 'Ret - Ret';
save_name = 'Crosscompare_RSA_Matfiles/CrossComp_RSA_Ret_Ret_btROI_noperm_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1:3
    jobHandles{end_ind+msr} = batch(c, @Crosscomparison_RSA_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
end

jobHandles{1}
diary(jobHandles{1})
%delete(jobHandles{1})


% Encoding to Retrieval, Within ROI, PermTest
cfg.permtest         = false; %true;
cfg.timewind2        = [-2.5 0.2];
cfg.ROI              = {'Wholehead'};
% cfg.ROI              = {'OCC';...
%                         'TMP';...
%                         'FRT';...
%                         'PRT'};

pair = 'Enc - Ret';
save_name = 'Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Ret_wiROI_noperm_wholehead_msr';
c = parcluster();
%jobHandles = {};
end_ind = length(jobHandles);
for msr = 1:3
    jobHandles{end_ind+msr} = batch(c, @Crosscomparison_RSA_Enc_Ret_LoadData, 0, {msr, cfg, pair, save_name},'Pool',5);
end

jobHandles{4}
diary(jobHandles{4})




%% Plot Crosscomparison with Encoding and Retrieval Time series

% Subject Names
load('RSA_Data_Enc','Subj_names')

% Create Data Struct for Encoding
RSA_Data_Enc = [];
RSA_Data_Ret = [];
for sub = 1:length(Subj_names)
    tmp_strct_enc = load(['RSA_Data/RSA_Data_Enc_',Subj_names{sub}]);
    tmp_strct_ret = load(['RSA_Data/RSA_Data_Ret_',Subj_names{sub}]);
    if(sub == 1)
        ROI = fieldnames(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]));
        RSA_Data_Enc.Names   = Subj_names;
        for fn = fieldnames(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{1}))'
            RSA_Data_Enc.(fn{1}) = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).(fn{1});
        end
        RSA_Data_Enc = rmfield(RSA_Data_Enc, {'Name','TimeVec1024','RSA_full','MDS_full','RSA_16','Encoding_Data','TrialInfo'...
                                      'rsa_dim','curROI','curROI_name'});
        for r = 1:length(ROI)
            RSA_Data_Enc.([ROI{r},'_ROI']) = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).curROI;
        end

        RSA_Data_Ret.Names   = Subj_names;
        for fn = fieldnames(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{1}))'
            RSA_Data_Ret.(fn{1}) = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).(fn{1});
        end
        RSA_Data_Ret = rmfield(RSA_Data_Ret, {'Name','TimeVec1024','RSA_full','MDS_full','RSA_16','Retrieval_Data','TrialInfo'...
                                      'rsa_dim','curROI','curROI_name'});
        for r = 1:length(ROI)
            RSA_Data_Ret.([ROI{r},'_ROI']) = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).curROI;
        end

    end
    
    for r = 1:length(ROI)
        if(~isempty(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).RSA_full))
            RSA_Data_Enc.(ROI{r}).full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_full{1,msr},[3 1 2]);
        end
        RSA_Data_Enc.(ROI{r}).red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_16{1,msr},[3 1 2]);
        
        if(~isempty(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).RSA_full))
            RSA_Data_Ret.(ROI{r}).full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_full{1,msr},[3 1 2]);
        end
        RSA_Data_Ret.(ROI{r}).red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_16{1,msr},[3 1 2]);
    end
    
    RSA_Data_Enc.Encoding_Data{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).Encoding_Data;
    RSA_Data_Enc.TrialInfo{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).TrialInfo;

    RSA_Data_Ret.Retrieval_Data{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).Retrieval_Data;
    RSA_Data_Ret.TrialInfo{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).TrialInfo;

end
clear tmp_strct_enc tmp_strct_ret



measures = {'LDA','SVM','euclidian','euclidian_wcc'};

tmp_struct = load('Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Ret_btROI_perm_msr1',['CrossComp_RSA_',measures{msr}]);
tmp_fnames = fieldnames(tmp_struct);

method = {'Corr','Meth1_per','Meth1_sem','Meth2_per','Meth2_sem'};
permtestname = {'CorrPermTest', 'Meth1_per_PermTest', 'Meth1_sem_PermTest',  'Meth2_per_PermTest', 'Meth2_sem_PermTest'};
mth1 = 2; mth2 = 3;

ROI = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).ROI;
if(size(ROI,2) == 2)
    ROI_names = strcat(ROI(:,1),'_',ROI(:,2));
else
    ROI_names = ROI;
end
Cat = {'Perceptual','Semantic'};
Col = {'b','r'};

time_wind_enc = [-0.1  1.2];
time_wind_ret = [-2.2 -0.8]; % [-1.8 -0.6]

r1 = 1; r2 = 2;
c1 = 1; c2 = 2;
dt1 = [3 4]; dt2 = [3 4];
rndcut_perm = true;

TimeX = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).TimeVec1;
time_idx_enc = dsearchn(TimeX',time_wind_enc')';
TimeY = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).TimeVec2;
time_idx_ret = dsearchn(TimeY',time_wind_ret')';
CrossComp_Data.RSA_red16.Data1 = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).RSA_red16.(ROI_names{r1}).(method{mth1})(:,time_idx_ret(1):time_idx_ret(end),time_idx_enc(1):time_idx_enc(end));
CrossComp_Data.RSA_red16.Data2 = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).RSA_red16.(ROI_names{r2}).(method{mth2})(:,time_idx_ret(1):time_idx_ret(end),time_idx_enc(1):time_idx_enc(end));
TimeX = TimeX(time_idx_enc(1):time_idx_enc(end));
TimeY = TimeY(time_idx_ret(1):time_idx_ret(end));


if(sum(dt1) == 7)
    Hyp_Mat{1} = double(Perceptual_Mat_red16>0)  + -double(Perceptual_Mat_red16<0);
    Hyp_Mat{2} = double(Semantic_Mat_red16>0)    + -double(Semantic_Mat_red16<0);
else
    Hyp_Mat{1} = double(Perceptual_Mat_red16==1) + -double(Perceptual_Mat_red16==2);
    Hyp_Mat{2} = double(Semantic_Mat_red16==1)   + -double(Semantic_Mat_red16==2);
end

cfg = [];
cfg.nPerms = 1000;
cfg.thresh_pval = 0.05;
cfg.mcc_cluster_pval = 0.05;
cfg.matshuff = false;
cfg.twoside = true;

cfg.Hyp_Mat = Hyp_Mat{c1};
cfg.TimeVec = RSA_Data_Enc.TimeVec(dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(1)):dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(end)));
Results1 = rsa_perm(cfg, RSA_Data_Enc.(ROI{r1,1}).red16_Data(:,dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(1)):dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(end)),:,:));
cfg.TimeVec = RSA_Data_Ret.TimeVec(dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(1)):dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(end)));
Results2 = rsa_perm(cfg, RSA_Data_Ret.(ROI{r1,end}).red16_Data(:,dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(1)):dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(end)),:,:));

cfg.Hyp_Mat = Hyp_Mat{c2};
cfg.TimeVec = RSA_Data_Enc.TimeVec(dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(1)):dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(end)));
Results3 = rsa_perm(cfg, RSA_Data_Enc.(ROI{r2,1}).red16_Data(:,dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(1)):dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(end)),:,:));
cfg.TimeVec = RSA_Data_Ret.TimeVec(dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(1)):dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(end)));
Results4 = rsa_perm(cfg, RSA_Data_Ret.(ROI{r2,end}).red16_Data(:,dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(1)):dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(end)),:,:));


if(rndcut_perm)
    cfg = [];
    cfg.nPerms = 1000;
    cfg.thresh_pval = 0.05;
    cfg.mcc_cluster_pval = 0.05;
    cfg.twoside = false;
    Results5 = PermTestRNDCP2D(cfg, CrossComp_Data.RSA_red16.Data1);
    Results6 = PermTestRNDCP2D(cfg, CrossComp_Data.RSA_red16.Data2);
else
    Results5.zmapthresh = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).(permtestname{mth1}).(ROI_names{r1}).zmapthresh;
    Results5.H = sum(~isnan(Results5.zmapthresh(:))) ~= 0;
    Results6.zmapthresh = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).(permtestname{mth2}).(ROI_names{r2}).zmapthresh;
    Results6.H = sum(~isnan(Results6.zmapthresh(:))) ~= 0;
end


enc_lim = [TimeX(1) TimeX(end)]; %enc_lim = [-0.2 1.0];
ret_lim = [TimeY(1) TimeY(end)]; %ret_lim = [-2.5 0.2];

figure('Pos',[189 137 1531 782])

h1 = subplot(1,2,1);
contourf(TimeX, TimeY, squeeze(nanmean(CrossComp_Data.RSA_red16.Data1,1)), 40,'linestyle','none'); colorbar
caxis([-.0055 .0055]); set(h1,'xlim', enc_lim, 'ylim', ret_lim, 'xticklabel',[], 'yticklabel',[]); 
title(ROI_names{r1});
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
if(Results5.H == 1)
    zmapthresh = Results5.zmapthresh; zmapthresh(isnan(zmapthresh)) = 0; zmapthresh(zmapthresh ~= 0) = 1;
    contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
end
hold off

aspect = get(h1,'PlotBoxAspectRatio');
set(h1,'Units','pixels');
pos = get(h1,'Position')+[0 100 0 -100];
pos(3) = aspect(1)/aspect(2)*pos(4);
set(h1,'Position',pos);

h2 = axes('Units','pixels','pos',[pos(1) pos(2)-130 pos(3) 128]);%,'visible','off');
Dat_names1 = fieldnames(RSA_Time_Enc.(ROI{r1,1}).(Cat{c1}));
dat1 = RSA_Time_Enc.(ROI{r1,1}).(Cat{c1}).(Dat_names1{dt1(2)})(:,dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(1)):dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(end))) - RSA_Time_Enc.(ROI{r1,1}).(Cat{c1}).(Dat_names1{dt1(1)})(:,dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(1)):dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(end)));
SEM1 = nanstd(dat1,0,1)./sqrt(size(dat1,1)); hold on;
fill([Results1.TimeVec fliplr(Results1.TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) + SEM1)],Col{c1},'FaceAlpha',0.3,'EdgeAlpha',0);
fill([Results1.TimeVec fliplr(Results1.TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) - SEM1)],Col{c1},'FaceAlpha',0.3,'EdgeAlpha',0);
plot([0 0],[min(nanmean(dat1,1) - SEM1) max(nanmean(dat1,1) + SEM1)]*1.4,'--k','linewidth',2)
l1 = plot(Results1.TimeVec, nanmean(dat1,1),Col{c1},'linewidth',2);
sign_mcc_clust_1 = Results1.zmapthresh;
sign_mcc_clust_1(abs(sign_mcc_clust_1) > 0) = min(nanmean(dat1,1) - SEM1);
plot(Results1.TimeVec,sign_mcc_clust_1,[Col{c1},'o'],'MarkerFaceColor',Col{c1})
ylabel(['Diff ',RSA_Data_Enc.meas16{msr}]); xlim(enc_lim); ylim([min(nanmean(dat1,1) - SEM1) max(nanmean(dat1,1) + SEM1)]*1.4); hold off
lg = legend(l1, sprintf('%s \n%s - %s',Cat{c1},Dat_names1{dt1(2)},Dat_names1{dt1(1)})); legend boxoff; set(lg,'FontSize',8)
%lg = legend(l1, {[Cat{c1},' \n ',Dat_names1{dt1(2)},' - ',Dat_names1{dt1(1)}]}); legend boxoff; set(lg,'FontSize',8)
xlabel('Encoding')

h3 = axes('Units','pixels','pos',[pos(1)-135 pos(2) 133 pos(4)]);%,'visible','off');
Dat_names1 = fieldnames(RSA_Time_Ret.(ROI{r1,end}).(Cat{c1}));
dat1 = RSA_Time_Ret.(ROI{r1,end}).(Cat{c1}).(Dat_names1{dt1(2)})(:,dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(1)):dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(end))) - RSA_Time_Ret.(ROI{r1,end}).(Cat{c1}).(Dat_names1{dt1(1)})(:,dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(1)):dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(end)));
SEM1 = nanstd(dat1,0,1)./sqrt(size(dat1,1)); hold on;
fill([nanmean(dat1,1) fliplr(nanmean(dat1,1) + SEM1)],[Results2.TimeVec fliplr(Results2.TimeVec)],Col{c1},'FaceAlpha',0.3,'EdgeAlpha',0);
fill([nanmean(dat1,1) fliplr(nanmean(dat1,1) - SEM1)],[Results2.TimeVec fliplr(Results2.TimeVec)],Col{c1},'FaceAlpha',0.3,'EdgeAlpha',0);
plot([min(nanmean(dat1,1) - SEM1) max(nanmean(dat1,1) + SEM1)]*1.4,[0 0],'--k','linewidth',2)
plot(nanmean(dat1,1), Results2.TimeVec,Col{c1},'linewidth',2);
sign_mcc_clust_2 = Results2.zmapthresh;
sign_mcc_clust_2(abs(sign_mcc_clust_2) > 0) = min(nanmean(dat1,1) - SEM1);
plot(sign_mcc_clust_2,Results2.TimeVec,[Col{c1},'o'],'MarkerFaceColor',Col{c1})
xlabel(['Diff ',RSA_Data_Ret.meas16{msr}]); ylim(ret_lim); xlim([min(nanmean(dat1,1) - SEM1) max(nanmean(dat1,1) + SEM1)]*1.4); hold off
set(h3,'xdir','reverse')
ylabel('Retrieval')


h4 = subplot(1,2,2);
contourf(TimeX, TimeY, squeeze(nanmean(CrossComp_Data.RSA_red16.Data2,1)), 40,'linestyle','none'); colorbar
caxis([-.0055 .0055]); set(h4,'xlim', enc_lim, 'ylim', ret_lim, 'xticklabel',[], 'yticklabel',[]); 
title(ROI_names{r2});
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
if(Results6.H == 1)
    zmapthresh = Results6.zmapthresh; zmapthresh(isnan(zmapthresh)) = 0; zmapthresh(zmapthresh ~= 0) = 1;
    contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
end
hold off

aspect = get(h4,'PlotBoxAspectRatio');
set(h4,'Units','pixels');
pos = get(h4,'Position')+[0 100 0 -100];
pos(3) = aspect(1)/aspect(2)*pos(4);
set(h4,'Position',pos);

h5 = axes('Units','pixels','pos',[pos(1) pos(2)-130 pos(3) 128]);%,'visible','off');
Dat_names1 = fieldnames(RSA_Time_Enc.(ROI{r2,1}).(Cat{c2}));
dat1 = RSA_Time_Enc.(ROI{r2,1}).(Cat{c2}).(Dat_names1{dt1(2)})(:,dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(1)):dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(end))) - RSA_Time_Enc.(ROI{r2,1}).(Cat{c2}).(Dat_names1{dt1(1)})(:,dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(1)):dsearchn(RSA_Data_Enc.TimeVec', time_wind_enc(end)));
SEM1 = nanstd(dat1,0,1)./sqrt(size(dat1,1)); hold on;
fill([Results3.TimeVec fliplr(Results3.TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) + SEM1)],Col{c2},'FaceAlpha',0.3,'EdgeAlpha',0);
fill([Results3.TimeVec fliplr(Results3.TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) - SEM1)],Col{c2},'FaceAlpha',0.3,'EdgeAlpha',0);
plot([0 0],[min(nanmean(dat1,1) - SEM1) max(nanmean(dat1,1) + SEM1)]*1.4,'--k','linewidth',2)
l1 = plot(Results3.TimeVec, nanmean(dat1,1),Col{c2},'linewidth',2);
sign_mcc_clust_1 = Results3.zmapthresh;
sign_mcc_clust_1(abs(sign_mcc_clust_1) > 0) = min(nanmean(dat1,1) - SEM1);
plot(Results3.TimeVec,sign_mcc_clust_1,[Col{c2},'o'],'MarkerFaceColor',Col{c2})
ylabel(['Diff ',RSA_Data_Enc.meas16{msr}]); xlim(enc_lim); ylim([min(nanmean(dat1,1) - SEM1) max(nanmean(dat1,1) + SEM1)]*1.4); hold off
lg = legend(l1, sprintf('%s \n%s - %s',Cat{c2},Dat_names1{dt1(2)},Dat_names1{dt1(1)})); legend boxoff; set(lg,'FontSize',8)
%lg = legend(l1, {[Cat{c1},' \n ',Dat_names1{dt1(2)},' - ',Dat_names1{dt1(1)}]}); legend boxoff; set(lg,'FontSize',8)
xlabel('Encoding')

h6 = axes('Units','pixels','pos',[pos(1)-135 pos(2) 133 pos(4)]);%,'visible','off');
Dat_names1 = fieldnames(RSA_Time_Ret.(ROI{r2,end}).(Cat{c2}));
dat1 = RSA_Time_Ret.(ROI{r2,end}).(Cat{c2}).(Dat_names1{dt1(2)})(:,dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(1)):dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(end))) - RSA_Time_Ret.(ROI{r2,end}).(Cat{c2}).(Dat_names1{dt1(1)})(:,dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(1)):dsearchn(RSA_Data_Ret.TimeVec', time_wind_ret(end)));
SEM1 = nanstd(dat1,0,1)./sqrt(size(dat1,1)); hold on;
fill([nanmean(dat1,1) fliplr(nanmean(dat1,1) + SEM1)],[Results4.TimeVec fliplr(Results4.TimeVec)],Col{c2},'FaceAlpha',0.3,'EdgeAlpha',0);
fill([nanmean(dat1,1) fliplr(nanmean(dat1,1) - SEM1)],[Results4.TimeVec fliplr(Results4.TimeVec)],Col{c2},'FaceAlpha',0.3,'EdgeAlpha',0);
plot([min(nanmean(dat1,1) - SEM1) max(nanmean(dat1,1) + SEM1)]*1.4,[0 0],'--k','linewidth',2)
plot(nanmean(dat1,1), Results4.TimeVec,Col{c2},'linewidth',2);
sign_mcc_clust_2 = Results4.zmapthresh;
sign_mcc_clust_2(abs(sign_mcc_clust_2) > 0) = min(nanmean(dat1,1) - SEM1);
plot(sign_mcc_clust_2,Results4.TimeVec,[Col{c2},'o'],'MarkerFaceColor',Col{c2})
xlabel(['Diff ',RSA_Data_Ret.meas16{msr}]); ylim(ret_lim); xlim([min(nanmean(dat1,1) - SEM1) max(nanmean(dat1,1) + SEM1)]*1.4); hold off
set(h6,'xdir','reverse')
ylabel('Retrieval')

set([h1 h2 h3 h4 h5 h6],'Units','normalized')




%% Plot all ROIs and Methods


measures = {'LDA','SVM','euclidian','euclidian_wcc'};
msr = 1;
method = {'Corr','Meth1_per','Meth1_sem','Meth2_per','Meth2_sem'};
mth1 = 2; mth2 = 3;



%% Encoding to Encoding  Within ROIs

tmp_struct = load(['Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Enc_wiROI_noperm_msr',num2str(msr)],['CrossComp_RSA_',measures{msr}]);
Data = tmp_struct.(['CrossComp_RSA_',measures{msr}]).RSA_red16;
TimeX = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec1;
TimeY = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec2;
ROI = tmp_struct.(['CrossComp_RSA_',measures{msr}]).ROI;

enc_lim = [-0.2  1.5];


figure('Pos',[484 44 885 952])
ct = 1;
for m = 1:length(method)
    for r = 1:size(ROI,1)
        curData = squeeze(nanmean(Data.(ROI{r,1}).(method{m}),1)); prctile(curData(:),[1 99]);
        subplot(length(method),size(ROI,1),ct)
        contourf(TimeX, TimeY, curData, 40,'linestyle','none'); ax_pos = get(gca,'Pos');
        cl = colorbar('Position',[ax_pos(1)+ax_pos(3)+0.005 ax_pos(2) 0.01 ax_pos(4)]);
        caxis([min(curData(:)) prctile(curData(:),80)]); set(gca,'xlim', enc_lim, 'ylim', enc_lim, 'xticklabel',[], 'yticklabel',[]); 
        hold on
        plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
        plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
        hold off
        ct = ct + 1;
        if(m == 1)
            text('string',ROI{r,1},'Pos',[0.5*(enc_lim(end) + enc_lim(1)) enc_lim(end)*1.2 0],'FontSize',11,'HorizontalAlignment','center')
        end
        if(r == 1)
            text('string',strrep(method{m},'_',' '),'Pos',[enc_lim(1)*2 0.5*(enc_lim(end) + enc_lim(1)) 0],'FontSize',11,'HorizontalAlignment','center','Rotation',90)
        end
    end
end


r = 4; m = 1;
figure
curData = squeeze(nanmean(Data.(ROI{r}).(method{m}),1)); prctile(curData(:),[1 99]);
contourf(TimeX, TimeY, curData, 40,'linestyle','none'); cl = colorbar;
caxis([min(curData(:)) prctile(curData(:),70)]); set(gca,'xlim', enc_lim, 'ylim', enc_lim);
xlabel('Encoding'); ylabel('Encoding'); title(ROI{r}); ylabel(cl, strrep(method{m},'_',' '))
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
hold off
saveas(gcf,sprintf('Results/Enc_Enc_%s_%s_TGM_%s.png',ROI{r},method{m},measures{msr}))


%% Encoding to Encoding  Between ROIs

tmp_struct = load(['Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Enc_btROI_noperm_msr',num2str(msr)],['CrossComp_RSA_',measures{msr}]);
Data = tmp_struct.(['CrossComp_RSA_',measures{msr}]).RSA_red16;
TimeX = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec1;
TimeY = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec2;
ROI = tmp_struct.(['CrossComp_RSA_',measures{msr}]).ROI;

enc1_lim = [-0.2  1];
enc2_lim = [-0.2  1];


figure('Pos',[300 42 1267 952])
ct = 1;
for m = 1:length(method)
    for r = 1:size(ROI,1)
        curData = squeeze(nanmean(Data.([ROI{r,1},'_',ROI{r,2}]).(method{m}),1)); prctile(curData(:),[1 99]);
        subplot(length(method),size(ROI,1),ct)
        contourf(TimeX, TimeY, curData, 40,'linestyle','none'); ax_pos = get(gca,'Pos');
        cl = colorbar('Position',[ax_pos(1)+ax_pos(3)+0.005 ax_pos(2) 0.01 ax_pos(4)]);
        caxis(prctile(curData(:),[2.5 97.5])); set(gca,'xlim', enc1_lim, 'ylim', enc2_lim, 'xticklabel',[], 'yticklabel',[]); 
        hold on
        plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
        plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
        hold off
        ct = ct + 1;
        if(m == 1)
            text('string',[ROI{r,1},' - ',ROI{r,2}],'Pos',[0.5*(enc1_lim(end) + enc1_lim(1)) enc2_lim(end)*1.2 0],'FontSize',11,'HorizontalAlignment','center')
        end
        if(r == 1)
            text('string',strrep(method{m},'_',' '),'Pos',[enc1_lim(1)*2 0.5*(enc2_lim(end) + enc2_lim(1)) 0],'FontSize',11,'HorizontalAlignment','center','Rotation',90)
        end
    end
end


r = 1; m = 1;
figure
curData = squeeze(nanmean(Data.([ROI{r,1},'_',ROI{r,2}]).(method{m}),1)); prctile(curData(:),[1 99]);
contourf(TimeX, TimeY, curData, 40,'linestyle','none'); cl = colorbar;
caxis(prctile(curData(:),[2.5 97.5])); set(gca,'xlim', enc1_lim, 'ylim', enc2_lim); 
xlabel(['Encoding ',ROI{r,1}]); ylabel(['Encoding ',ROI{r,2}]); ylabel(cl, strrep(method{m},'_',' '))
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
hold off



%% Retrieval to Retrieval  Within ROIs

tmp_struct = load(['Crosscompare_RSA_Matfiles/CrossComp_RSA_Ret_Ret_wiROI_noperm_msr',num2str(msr)],['CrossComp_RSA_',measures{msr}]);
Data = tmp_struct.(['CrossComp_RSA_',measures{msr}]).RSA_red16;
TimeX = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec1;
TimeY = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec2;
ROI = tmp_struct.(['CrossComp_RSA_',measures{msr}]).ROI;

ret_lim = [-1.8  -0.6];


figure('Pos',[484 44 885 952])
ct = 1;
for m = 1:length(method)
    for r = 1:size(ROI,1)
        curData = squeeze(nanmean(Data.(ROI{r,1}).(method{m}),1)); prctile(curData(:),[1 99]);
        subplot(length(method),size(ROI,1),ct)
        contourf(TimeX, TimeY, curData, 40,'linestyle','none'); ax_pos = get(gca,'Pos');
        cl = colorbar('Position',[ax_pos(1)+ax_pos(3)+0.005 ax_pos(2) 0.01 ax_pos(4)]);
        caxis([min(curData(:)) prctile(curData(:),90)]); set(gca,'xlim', ret_lim, 'ylim', ret_lim, 'xticklabel',[], 'yticklabel',[]); 
        hold on
        plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
        plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
        hold off
        ct = ct + 1;
        if(m == 1)
            text('string',ROI{r,1},'Pos',[0.5*(ret_lim(end) + ret_lim(1)) ret_lim(end)*3 0],'FontSize',11,'HorizontalAlignment','center')
        end
        if(r == 1)
            text('string',strrep(method{m},'_',' '),'Pos',[ret_lim(1)*1.1 0.5*(ret_lim(end) + ret_lim(1)) 0],'FontSize',11,'HorizontalAlignment','center','Rotation',90)
        end
    end
end


r = 1; m = 1;
figure
curData = squeeze(nanmean(Data.(ROI{r}).(method{m}),1)); prctile(curData(:),[1 99]);
contourf(TimeX, TimeY, curData, 40,'linestyle','none'); cl = colorbar;
caxis([min(curData(:)) prctile(curData(:),90)]); set(gca,'xlim', ret_lim, 'ylim', ret_lim);
xlabel('Retrieval'); ylabel('Retrieval'); title(ROI{r}); ylabel(cl, strrep(method{m},'_',' '))
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
hold off



%% Retrieval to Retrieval  Between ROIs

tmp_struct = load(['Crosscompare_RSA_Matfiles/CrossComp_RSA_Ret_Ret_btROI_noperm_msr',num2str(msr)],['CrossComp_RSA_',measures{msr}]);
Data = tmp_struct.(['CrossComp_RSA_',measures{msr}]).RSA_red16;
TimeX = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec1;
TimeY = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec2;
ROI = tmp_struct.(['CrossComp_RSA_',measures{msr}]).ROI;

ret1_lim = [-1.8 -0.6];
ret2_lim = [-1.8 -0.6];


figure('Pos',[484 44 885 952])
ct = 1;
for m = 1:length(method)
    for r = 1:size(ROI,1)
        curData = squeeze(nanmean(Data.([ROI{r,1},'_',ROI{r,2}]).(method{m}),1)); prctile(curData(:),[1 99]);
        subplot(length(method),size(ROI,1),ct)
        contourf(TimeX, TimeY, curData, 40,'linestyle','none'); ax_pos = get(gca,'Pos');
        cl = colorbar('Position',[ax_pos(1)+ax_pos(3)+0.005 ax_pos(2) 0.01 ax_pos(4)]);
        caxis(prctile(curData(:),[2.5 97.5])); set(gca,'xlim', ret1_lim, 'ylim', ret2_lim, 'xticklabel',[], 'yticklabel',[]); 
        hold on
        plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
        plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
        hold off
        ct = ct + 1;
        if(m == 1)
            text('string',[ROI{r,1},' - ',ROI{r,2}],'Pos',[0.5*(ret1_lim(end) + ret1_lim(1)) ret2_lim(end)*1.2 0],'FontSize',11,'HorizontalAlignment','center')
        end
        if(r == 1)
            text('string',strrep(method{m},'_',' '),'Pos',[ret1_lim(1)*1.2 0.5*(ret2_lim(end) + ret2_lim(1)) 0],'FontSize',11,'HorizontalAlignment','center','Rotation',90)
        end
    end
end


r = 6; m = 2;
figure
curData = squeeze(nanmean(Data.([ROI{r,1},'_',ROI{r,2}]).(method{m}),1)); prctile(curData(:),[1 99]);
contourf(TimeX, TimeY, curData, 40,'linestyle','none'); cl = colorbar;
caxis(prctile(curData(:),[2.5 97.5])); set(gca,'xlim', ret1_lim, 'ylim', ret2_lim); 
xlabel(['Retrieval ',ROI{r,1}]); ylabel(['Retrieval ',ROI{r,2}]); ylabel(cl, strrep(method{m},'_',' '))
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
hold off



%% Encoding to Retrieval  Within ROIs

tmp_struct = load('Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Ret_wiROI_noperm_msr1',['CrossComp_RSA_',measures{msr}]);
Data = tmp_struct.(['CrossComp_RSA_',measures{msr}]).RSA_red16;
TimeX = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec1;
TimeY = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec2;
ROI = tmp_struct.(['CrossComp_RSA_',measures{msr}]).ROI;

enc_lim = [-0.1  1.5];
ret_lim = [-1.5  -0.5];


figure('Pos',[484 44 885 952])
ct = 1;
for m = 1:length(method)
    for r = 1:size(ROI,1)
        curData = squeeze(nanmean(Data.(ROI{r,1}).(method{m}),1)); prctile(curData(:),[1 99]);
        subplot(length(method),size(ROI,1),ct)
        contourf(TimeX, TimeY, curData, 40,'linestyle','none'); %colorbar
        caxis(prctile(curData(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim, 'xticklabel',[], 'yticklabel',[]); 
        hold on
        plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
        plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
        hold off
        ct = ct + 1;
        if(m == 1)
            text('string',ROI{r,1},'Pos',[0.5*(enc_lim(end) + enc_lim(1)) ret_lim(end)*1.2 0],'FontSize',11,'HorizontalAlignment','center')
        end
        if(r == 1)
            text('string',strrep(method{m},'_',' '),'Pos',[enc_lim(1)*2 0.5*(ret_lim(end) + ret_lim(1)) 0],'FontSize',11,'HorizontalAlignment','center','Rotation',90)
        end
    end
end


r = 4; m = 3;
figure
curData = squeeze(nanmean(Data.(ROI{r}).(method{m}),1)); prctile(curData(:),[1 99]);
contourf(TimeX, TimeY, curData, 40,'linestyle','none'); cl = colorbar;
caxis(prctile(curData(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim);
xlabel('Encoding'); ylabel('Retrieval'); title(ROI{r}); ylabel(cl, strrep(method{m},'_',' '))
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
hold off



%% Encoding to Retrieval  Between ROIs

tmp_struct = load(['Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Ret_btROI_noperm_msr',num2str(msr)],['CrossComp_RSA_',measures{msr}]);
Data = tmp_struct.(['CrossComp_RSA_',measures{msr}]).RSA_red16;
TimeX = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec1;
TimeY = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec2;
ROI = tmp_struct.(['CrossComp_RSA_',measures{msr}]).ROI;

enc_lim = [-0.1  1.2];
ret_lim = [-2.4 -1.6];


figure('Pos',[484 44 885 952])
ct = 1;
for m = 1:length(method)-2
    for r = 1:size(ROI,1)
        curData = squeeze(nanmean(Data.([ROI{r,1},'_',ROI{r,2}]).(method{m}),1)); prctile(curData(:),[1 99]);
        subplot(length(method)-2,size(ROI,1),ct)
        contourf(TimeX, TimeY, curData, 40,'linestyle','none'); %colorbar
        caxis(prctile(curData(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim, 'xticklabel',[], 'yticklabel',[]); 
        hold on
        plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
        plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
        hold off
        ct = ct + 1;
        if(m == 1)
            text('string',[ROI{r,1},' - ',ROI{r,2}],'Pos',[0.5*(enc_lim(end) + enc_lim(1)) ret_lim(end)*1.2 0],'FontSize',11,'HorizontalAlignment','center')
        end
        if(r == 1)
            text('string',strrep(method{m},'_',' '),'Pos',[enc_lim(1)*2 0.5*(ret_lim(end) + ret_lim(1)) 0],'FontSize',11,'HorizontalAlignment','center','Rotation',90)
        end
    end
end


r = 3; m = 2;
figure
curData = squeeze(nanmean(Data.([ROI{r,1},'_',ROI{r,2}]).(method{m}),1)); prctile(curData(:),[1 99]);
contourf(TimeX, TimeY, curData, 40,'linestyle','none'); cl = colorbar;
caxis(prctile(curData(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim); 
xlabel(['Encoding ',ROI{r,1}]); ylabel(['Retrieval ',ROI{r,2}]); ylabel(cl, strrep(method{m},'_',' '))
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
hold off




%% Encoding to Retrieval  Within ROIs with ROIs

measures = {'LDA','SVM','euclidian','euclidian_wcc'};

tmp_struct = load('Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Ret_wiROI_noperm_msr',['CrossComp_RSA_',measures{msr}]);
Data = tmp_struct.(['CrossComp_RSA_',measures{msr}]).RSA_red16;
TimeX = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec1;
TimeY = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec2;
ROI = tmp_struct.(['CrossComp_RSA_',measures{msr}]).ROI;

enc_lim = [-0.1  1.2];
ret_lim = [-1.8  -0.6];

tmp_fnames = fieldnames(tmp_struct);

method = {'Corr','Meth1_per','Meth1_sem','Meth2_per','Meth2_sem'};
permtestname = {'CorrPermTest', 'Meth1_per_PermTest', 'Meth1_sem_PermTest',  'Meth2_per_PermTest', 'Meth2_sem_PermTest'};
mth1 = 2; mth2 = 3;

ROI = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).ROI;
if(size(ROI,2) == 2)
    ROI_names = strcat(ROI(:,1),'_',ROI(:,2));
else
    ROI_names = ROI;
end
Cat = {'Perceptual','Semantic'};
Col = {'b','r'};

r1 = 4; r2 = 4;
c1 = 1; c2 = 2;
cp = 1;


Results5.zmapthresh = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).(permtestname{mth1}).(ROI_names{r1}).zmapthresh;
Results5.H = sum(~isnan(Results5.zmapthresh(:))) ~= 0;
Results6.zmapthresh = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).(permtestname{mth2}).(ROI_names{r2}).zmapthresh;
Results6.H = sum(~isnan(Results6.zmapthresh(:))) ~= 0;


ROI1 = {'Wholehead'}; r = 1;

elec_idx = ROI_idx{r};
Sensor_Data_Enc = zeros(size(RSA_Data_Enc.(ROI1{r}).red16_SensorData,1), length(ROI_all_idx), size(RSA_Data_Enc.(ROI1{r}).red16_SensorData,3),4);
Sensor_Data_Enc(:,elec_idx,:,1) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,1),4);
Sensor_Data_Enc(:,elec_idx,:,2) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,2),4);
Sensor_Data_Enc(:,elec_idx,:,3) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,3),4);
Sensor_Data_Enc(:,elec_idx,:,4) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,4),4);

GA_Sensor_Data_Enc = [];
GA_Sensor_Data_Enc.avg1_1 = squeeze(mean(Sensor_Data_Enc(:,:,:,1),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc.var1_1 = squeeze(std(Sensor_Data_Enc(:,:,:,1),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Enc.avg1_2 = squeeze(mean(Sensor_Data_Enc(:,:,:,2),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc.var1_2 = squeeze(std(Sensor_Data_Enc(:,:,:,2),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1)); 
GA_Sensor_Data_Enc.avg2_1 = squeeze(mean(Sensor_Data_Enc(:,:,:,3),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc.var2_1 = squeeze(std(Sensor_Data_Enc(:,:,:,3),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Enc.avg2_2 = squeeze(mean(Sensor_Data_Enc(:,:,:,4),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc.var2_2 = squeeze(std(Sensor_Data_Enc(:,:,:,4),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Enc.dof = size(Sensor_Data_Enc,1).*ones(size(Sensor_Data_Enc,2),size(Sensor_Data_Enc,3)); 
GA_Sensor_Data_Enc.time = RSA_Data_Enc.TimeVec; 
GA_Sensor_Data_Enc.label = elecs; 
GA_Sensor_Data_Enc.dimord = 'chan_time';

elec_idx = ROI_idx{r};
Sensor_Data_Ret = zeros(size(RSA_Data_Ret.(ROI1{r}).red16_SensorData,1), length(ROI_all_idx), size(RSA_Data_Ret.(ROI1{r}).red16_SensorData,3),4);
Sensor_Data_Ret(:,elec_idx,:,1) = mean(RSA_Data_Ret.(ROI1{r}).red16_SensorData(:,:,:,1),4);
Sensor_Data_Ret(:,elec_idx,:,2) = mean(RSA_Data_Ret.(ROI1{r}).red16_SensorData(:,:,:,2),4);
Sensor_Data_Ret(:,elec_idx,:,3) = mean(RSA_Data_Ret.(ROI1{r}).red16_SensorData(:,:,:,3),4);
Sensor_Data_Ret(:,elec_idx,:,4) = mean(RSA_Data_Ret.(ROI1{r}).red16_SensorData(:,:,:,4),4);

GA_Sensor_Data_Ret = [];
GA_Sensor_Data_Ret.avg1_1 = squeeze(mean(Sensor_Data_Ret(:,:,:,1),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Ret.var1_1 = squeeze(std(Sensor_Data_Ret(:,:,:,1),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Ret.avg1_2 = squeeze(mean(Sensor_Data_Ret(:,:,:,2),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Ret.var1_2 = squeeze(std(Sensor_Data_Ret(:,:,:,2),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1)); 
GA_Sensor_Data_Ret.avg2_1 = squeeze(mean(Sensor_Data_Ret(:,:,:,3),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Ret.var2_1 = squeeze(std(Sensor_Data_Ret(:,:,:,3),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Ret.avg2_2 = squeeze(mean(Sensor_Data_Ret(:,:,:,4),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Ret.var2_2 = squeeze(std(Sensor_Data_Ret(:,:,:,4),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Ret.dof = size(Sensor_Data_Ret,1).*ones(size(Sensor_Data_Ret,2),size(Sensor_Data_Ret,3)); 
GA_Sensor_Data_Ret.time = RSA_Data_Ret.TimeVec; 
GA_Sensor_Data_Ret.label = elecs; 
GA_Sensor_Data_Ret.dimord = 'chan_time';


cfg = [];
cfg.parameter        = [];
cfg.xlim             = []; %[0.24:0.02:0.54]; %[0.06:0.02:0.36];
cfg.zlim             = [];
%cfg.colorbar         = 'yes';
cfg.layout           = 'biosemi128.lay';
cfg.highlight        = 'on';
cfg.highlightchannel = elec_idx;
cfg.highlightsymbol  = 's';
cfg.highlightcolor   = [0 0 0];
cfg.comment          = 'no';


figure('Pos',[298  70  1285  924])
subplot(4,4,[1 2 5 6])
curData = squeeze(nanmean(Data.(ROI{r1}).(method{mth1}),1)); prctile(curData(:),[1 99]);
contourf(TimeX, TimeY, curData, 40,'linestyle','none'); cl = colorbar;
caxis(prctile(curData(:),[2.5 97.5])); set(gca,'xlim', enc_lim, 'ylim', ret_lim);
ylabel('Retrieval'); title(ROI{r1}); ylabel(cl, strrep(method{mth1},'_',' '))
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
if(Results5.H == 1)
    zmapthresh = Results5.zmapthresh; zmapthresh(isnan(zmapthresh)) = 0; zmapthresh(zmapthresh ~= 0) = 1;
    contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
end
hold off

subplot(4,4,[9 10 13 14])
curData = squeeze(nanmean(Data.(ROI{r2}).(method{mth2}),1)); prctile(curData(:),[1 99]);
contourf(TimeX, TimeY, curData, 40,'linestyle','none'); cl = colorbar;
caxis(prctile(curData(:),[2.5 97.5])); set(gca,'xlim', enc_lim, 'ylim', ret_lim);
xlabel('Encoding'); ylabel('Retrieval'); title(ROI{r2}); ylabel(cl, strrep(method{mth2},'_',' '))
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
if(Results6.H == 1)
    zmapthresh = Results6.zmapthresh; zmapthresh(isnan(zmapthresh)) = 0; zmapthresh(zmapthresh ~= 0) = 1;
    contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
end
hold off

subplot(4,4,3)
cfg.parameter = ['avg',num2str(c1),'_2'];
cfg.xlim      = [0.08 0.2];
cfg.zlim      = prctile(GA_Sensor_Data_Enc.(['avg',num2str(c1),'_2'])(:),[2.5 97.5]);
ft_topoplotER(cfg, GA_Sensor_Data_Enc);

subplot(4,4,4)
cfg.parameter = ['avg',num2str(c1),'_1'];
cfg.xlim      = [0.08 0.2];
cfg.zlim      = prctile(GA_Sensor_Data_Enc.(['avg',num2str(c1),'_1'])(:),[2.5 97.5]);
ft_topoplotER(cfg, GA_Sensor_Data_Enc);

subplot(4,4,7)
cfg.parameter = ['avg',num2str(c1),'_2'];
cfg.xlim      = [-1.3 -0.95];
cfg.zlim      = prctile(GA_Sensor_Data_Ret.(['avg',num2str(c1),'_2'])(:),[2.5 97.5]);
ft_topoplotER(cfg, GA_Sensor_Data_Ret);

subplot(4,4,8)
cfg.parameter = ['avg',num2str(c1),'_1'];
cfg.xlim      = [-1.3 -0.95];
cfg.zlim      = prctile(GA_Sensor_Data_Ret.(['avg',num2str(c1),'_1'])(:),[2.5 97.5]);
ft_topoplotER(cfg, GA_Sensor_Data_Ret);



subplot(4,4,11)
cfg.parameter = ['avg',num2str(c1),'_2'];
cfg.xlim      = [0.75 1.2];
cfg.zlim      = prctile(GA_Sensor_Data_Enc.(['avg',num2str(c1),'_2'])(:),[2.5 97.5]);
ft_topoplotER(cfg, GA_Sensor_Data_Enc);

subplot(4,4,12)
cfg.parameter = ['avg',num2str(c1),'_1'];
cfg.xlim      = [0.75 1.2];
cfg.zlim      = prctile(GA_Sensor_Data_Enc.(['avg',num2str(c1),'_1'])(:),[2.5 97.5]);
ft_topoplotER(cfg, GA_Sensor_Data_Enc);

subplot(4,4,15)
cfg.parameter = ['avg',num2str(c1),'_2'];
cfg.xlim      = [-1.45 -1.15];
cfg.zlim      = prctile(GA_Sensor_Data_Ret.(['avg',num2str(c1),'_2'])(:),[2.5 97.5]);
ft_topoplotER(cfg, GA_Sensor_Data_Ret);

subplot(4,4,16)
cfg.parameter = ['avg',num2str(c1),'_1'];
cfg.xlim      = [-1.45 -1.15];
cfg.zlim      = prctile(GA_Sensor_Data_Ret.(['avg',num2str(c1),'_1'])(:),[2.5 97.5]);
ft_topoplotER(cfg, GA_Sensor_Data_Ret);





%% Encoding to Retrieval  Between ROIs with ROIs

tmp_struct = load('Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Ret_btROI_perm_msr1',['CrossComp_RSA_',measures{msr}]);
Data = tmp_struct.(['CrossComp_RSA_',measures{msr}]).RSA_red16;
TimeX = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec1;
TimeY = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec2;
ROI = tmp_struct.(['CrossComp_RSA_',measures{msr}]).ROI;

enc_lim = [-0.2  1.5];
ret_lim = [-2.2  -0.8];

tmp_fnames = fieldnames(tmp_struct);

method = {'Corr','Meth1_per','Meth1_sem','Meth2_per','Meth2_sem'};
permtestname = {'CorrPermTest', 'Meth1_per_PermTest', 'Meth1_sem_PermTest',  'Meth2_per_PermTest', 'Meth2_sem_PermTest'};
mth1 = 2; mth2 = 3;

ROI = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).ROI;
if(size(ROI,2) == 2)
    ROI_names = strcat(ROI(:,1),'_',ROI(:,2));
else
    ROI_names = ROI;
end
Cat = {'Perceptual','Semantic'};
Col = {'b','r'};

r1 = 1; r2 = 2;
c1 = 1; c2 = 2;
cp = 1;


Results5.zmapthresh = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).(permtestname{mth1}).(ROI_names{r1}).zmapthresh;
Results5.H = sum(~isnan(Results5.zmapthresh(:))) ~= 0;
Results6.zmapthresh = tmp_struct.(tmp_fnames{~cellfun(@isempty, strfind(tmp_fnames,['CrossComp_RSA_',measures{msr}]))}).(permtestname{mth2}).(ROI_names{r2}).zmapthresh;
Results6.H = sum(~isnan(Results6.zmapthresh(:))) ~= 0;


ROI1 = {'Wholehead'}; r = 1;

elec_idx = ROI_idx{1}; %1
Sensor_Data_Enc1 = zeros(size(RSA_Data_Enc.(ROI1{r}).red16_SensorData,1), length(ROI_all_idx), size(RSA_Data_Enc.(ROI1{r}).red16_SensorData,3),4);
Sensor_Data_Enc1(:,elec_idx,:,1) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,1),4);
Sensor_Data_Enc1(:,elec_idx,:,2) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,2),4);
Sensor_Data_Enc1(:,elec_idx,:,3) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,3),4);
Sensor_Data_Enc1(:,elec_idx,:,4) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,4),4);

GA_Sensor_Data_Enc1 = [];
GA_Sensor_Data_Enc1.avg1_1 = squeeze(mean(Sensor_Data_Enc1(:,:,:,1),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc1.var1_1 = squeeze(std(Sensor_Data_Enc1(:,:,:,1),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Enc1.avg1_2 = squeeze(mean(Sensor_Data_Enc1(:,:,:,2),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc1.var1_2 = squeeze(std(Sensor_Data_Enc1(:,:,:,2),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1)); 
GA_Sensor_Data_Enc1.avg2_1 = squeeze(mean(Sensor_Data_Enc1(:,:,:,3),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc1.var2_1 = squeeze(std(Sensor_Data_Enc1(:,:,:,3),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Enc1.avg2_2 = squeeze(mean(Sensor_Data_Enc1(:,:,:,4),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc1.var2_2 = squeeze(std(Sensor_Data_Enc1(:,:,:,4),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Enc1.dof = size(Sensor_Data_Enc1,1).*ones(size(Sensor_Data_Enc1,2),size(Sensor_Data_Enc1,3)); 
GA_Sensor_Data_Enc1.time = RSA_Data_Enc.TimeVec; 
GA_Sensor_Data_Enc1.label = elecs; 
GA_Sensor_Data_Enc1.dimord = 'chan_time';


elec_idx = ROI_idx{r}; %2
Sensor_Data_Enc2 = zeros(size(RSA_Data_Enc.(ROI1{r}).red16_SensorData,1), length(ROI_all_idx), size(RSA_Data_Enc.(ROI1{r}).red16_SensorData,3),4);
Sensor_Data_Enc2(:,elec_idx,:,1) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,1),4);
Sensor_Data_Enc2(:,elec_idx,:,2) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,2),4);
Sensor_Data_Enc2(:,elec_idx,:,3) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,3),4);
Sensor_Data_Enc2(:,elec_idx,:,4) = mean(RSA_Data_Enc.(ROI1{r}).red16_SensorData(:,:,:,4),4);

GA_Sensor_Data_Enc2 = [];
GA_Sensor_Data_Enc2.avg1_1 = squeeze(mean(Sensor_Data_Enc2(:,:,:,1),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc2.var1_1 = squeeze(std(Sensor_Data_Enc2(:,:,:,1),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Enc2.avg1_2 = squeeze(mean(Sensor_Data_Enc2(:,:,:,2),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc2.var1_2 = squeeze(std(Sensor_Data_Enc2(:,:,:,2),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1)); 
GA_Sensor_Data_Enc2.avg2_1 = squeeze(mean(Sensor_Data_Enc2(:,:,:,3),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc2.var2_1 = squeeze(std(Sensor_Data_Enc2(:,:,:,3),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Enc2.avg2_2 = squeeze(mean(Sensor_Data_Enc2(:,:,:,4),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Enc2.var2_2 = squeeze(std(Sensor_Data_Enc2(:,:,:,4),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Enc2.dof = size(Sensor_Data_Enc2,1).*ones(size(Sensor_Data_Enc2,2),size(Sensor_Data_Enc2,3)); 
GA_Sensor_Data_Enc2.time = RSA_Data_Enc.TimeVec; 
GA_Sensor_Data_Enc2.label = elecs; 
GA_Sensor_Data_Enc2.dimord = 'chan_time';


elec_idx = ROI_idx{r}; %4
Sensor_Data_Ret = zeros(size(RSA_Data_Ret.(ROI1{r}).red16_SensorData,1), length(ROI_all_idx), size(RSA_Data_Ret.(ROI1{r}).red16_SensorData,3),4);
Sensor_Data_Ret(:,elec_idx,:,1) = mean(RSA_Data_Ret.(ROI1{r}).red16_SensorData(:,:,:,1),4);
Sensor_Data_Ret(:,elec_idx,:,2) = mean(RSA_Data_Ret.(ROI1{r}).red16_SensorData(:,:,:,2),4);
Sensor_Data_Ret(:,elec_idx,:,3) = mean(RSA_Data_Ret.(ROI1{r}).red16_SensorData(:,:,:,3),4);
Sensor_Data_Ret(:,elec_idx,:,4) = mean(RSA_Data_Ret.(ROI1{r}).red16_SensorData(:,:,:,4),4);

GA_Sensor_Data_Ret = [];
GA_Sensor_Data_Ret.avg1_1 = squeeze(mean(Sensor_Data_Ret(:,:,:,1),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Ret.var1_1 = squeeze(std(Sensor_Data_Ret(:,:,:,1),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Ret.avg1_2 = squeeze(mean(Sensor_Data_Ret(:,:,:,2),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Ret.var1_2 = squeeze(std(Sensor_Data_Ret(:,:,:,2),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1)); 
GA_Sensor_Data_Ret.avg2_1 = squeeze(mean(Sensor_Data_Ret(:,:,:,3),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Ret.var2_1 = squeeze(std(Sensor_Data_Ret(:,:,:,3),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Ret.avg2_2 = squeeze(mean(Sensor_Data_Ret(:,:,:,4),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data_Ret.var2_2 = squeeze(std(Sensor_Data_Ret(:,:,:,4),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data_Ret.dof = size(Sensor_Data_Ret,1).*ones(size(Sensor_Data_Ret,2),size(Sensor_Data_Ret,3)); 
GA_Sensor_Data_Ret.time = RSA_Data_Ret.TimeVec; 
GA_Sensor_Data_Ret.label = elecs; 
GA_Sensor_Data_Ret.dimord = 'chan_time';


cfg = [];
cfg.parameter        = [];
cfg.xlim             = []; %[0.24:0.02:0.54]; %[0.06:0.02:0.36];
cfg.zlim             = [];
%cfg.colorbar         = 'yes';
cfg.layout           = 'biosemi128.lay';
cfg.highlight        = 'on';
cfg.highlightchannel = elec_idx;
cfg.highlightsymbol  = 's';
cfg.highlightcolor   = [0 0 0];
cfg.comment          = 'no';


figure('Pos',[298  70  1285  924])
subplot(4,4,[1 2 5 6])
curData = squeeze(nanmean(Data.(ROI_names{r1}).(method{mth1}),1)); prctile(curData(:),[1 99]);
contourf(TimeX, TimeY, curData, 40,'linestyle','none'); cl = colorbar;
caxis(prctile(curData(:),[2.5 97.5])); set(gca,'xlim', enc_lim, 'ylim', ret_lim);
ylabel('Retrieval'); title(strrep(ROI_names{r1},'_',' ')); ylabel(cl, strrep(method{mth1},'_',' '))
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
if(Results5.H == 1)
    zmapthresh = Results5.zmapthresh; zmapthresh(isnan(zmapthresh)) = 0; zmapthresh(zmapthresh ~= 0) = 1;
    contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
end
hold off

subplot(4,4,[9 10 13 14])
curData = squeeze(nanmean(Data.(ROI_names{r2}).(method{mth2}),1)); prctile(curData(:),[1 99]);
contourf(TimeX, TimeY, curData, 40,'linestyle','none'); cl = colorbar;
caxis(prctile(curData(:),[2.5 97.5])); set(gca,'xlim', enc_lim, 'ylim', ret_lim);
xlabel('Encoding'); ylabel('Retrieval'); title(strrep(ROI_names{r2},'_',' ')); ylabel(cl, strrep(method{mth2},'_',' '))
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
if(Results6.H == 1)
    zmapthresh = Results6.zmapthresh; zmapthresh(isnan(zmapthresh)) = 0; zmapthresh(zmapthresh ~= 0) = 1;
    contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
end
hold off

subplot(4,4,3)
cfg.parameter = ['avg',num2str(c1),'_2'];
cfg.xlim      = [0.08 0.2];
cfg.zlim      = prctile(GA_Sensor_Data_Enc1.(['avg',num2str(c1),'_2'])(:),[2.5 97.5]);
cfg.highlightchannel = ROI_idx{1};
ft_topoplotER(cfg, GA_Sensor_Data_Enc1);

subplot(4,4,4)
cfg.parameter = ['avg',num2str(c1),'_1'];
cfg.xlim      = [0.08 0.2];
cfg.zlim      = prctile(GA_Sensor_Data_Enc1.(['avg',num2str(c1),'_1'])(:),[2.5 97.5]);
cfg.highlightchannel = ROI_idx{1};
ft_topoplotER(cfg, GA_Sensor_Data_Enc1);

subplot(4,4,7)
cfg.parameter = ['avg',num2str(c1),'_2'];
cfg.xlim      = [-1.3 -0.9];
cfg.zlim      = prctile(GA_Sensor_Data_Ret.(['avg',num2str(c1),'_2'])(:),[2.5 97.5]);
cfg.highlightchannel = ROI_idx{1}; %4
ft_topoplotER(cfg, GA_Sensor_Data_Ret);

subplot(4,4,8)
cfg.parameter = ['avg',num2str(c1),'_1'];
cfg.xlim      = [-1.3 -0.9];
cfg.zlim      = prctile(GA_Sensor_Data_Ret.(['avg',num2str(c1),'_1'])(:),[2.5 97.5]);
cfg.highlightchannel = ROI_idx{1}; %4
ft_topoplotER(cfg, GA_Sensor_Data_Ret);



subplot(4,4,11)
cfg.parameter = ['avg',num2str(c1),'_2'];
cfg.xlim      = [1.0 1.15];
cfg.zlim      = prctile(GA_Sensor_Data_Enc2.(['avg',num2str(c1),'_2'])(:),[2.5 97.5]);
cfg.highlightchannel = ROI_idx{1}; %2
ft_topoplotER(cfg, GA_Sensor_Data_Enc2);

subplot(4,4,12)
cfg.parameter = ['avg',num2str(c1),'_1'];
cfg.xlim      = [1.0 1.15];
cfg.zlim      = prctile(GA_Sensor_Data_Enc2.(['avg',num2str(c1),'_1'])(:),[2.5 97.5]);
cfg.highlightchannel = ROI_idx{1}; %2
ft_topoplotER(cfg, GA_Sensor_Data_Enc2);

subplot(4,4,15)
cfg.parameter = ['avg',num2str(c1),'_2'];
cfg.xlim      = [-2.15 -1.85];
cfg.zlim      = prctile(GA_Sensor_Data_Ret.(['avg',num2str(c1),'_2'])(:),[2.5 97.5]);
cfg.highlightchannel = ROI_idx{1}; %4
ft_topoplotER(cfg, GA_Sensor_Data_Ret);

subplot(4,4,16)
cfg.parameter = ['avg',num2str(c1),'_1'];
cfg.xlim      = [-2.15 -1.85];
cfg.zlim      = prctile(GA_Sensor_Data_Ret.(['avg',num2str(c1),'_1'])(:),[2.5 97.5]);
cfg.highlightchannel = ROI_idx{1}; %4
ft_topoplotER(cfg, GA_Sensor_Data_Ret);



