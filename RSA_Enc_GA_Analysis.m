%% RSA Encoding - All Subjects Analysis

%% Subject Names
load('RSA_Data_Enc','Subj_names')


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs

%% Measure
msr = 3;


%% Create Data
tmp_strct = load('RSA_Data_Enc');
RSA_Data = [];
for sub = 1:length(Subj_names)
    if(sub == 1)
        RSA_Data.Names   = Subj_names;
        for fn = fieldnames(tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC)'
            RSA_Data.(fn{1}) = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.(fn{1});
        end
        RSA_Data = rmfield(RSA_Data, {'Name','TimeVec1024','RSA_full','MDS_full','RSA_16','Encoding_Data','TrialInfo'...
                                      'rsa_dim','curROI','curROI_name'});
        RSA_Data.Encoding_Data{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.Encoding_Data;
        RSA_Data.TrialInfo{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;
        RSA_Data.OCC_ROI = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.curROI; 
        RSA_Data.TMP_ROI = tmp_strct.(['RSA_Data_',Subj_names{sub}]).TMP.curROI;
        RSA_Data.FRT_ROI = tmp_strct.(['RSA_Data_',Subj_names{sub}]).FRT.curROI;
        RSA_Data.CNT_ROI = tmp_strct.(['RSA_Data_',Subj_names{sub}]).CNT.curROI;
        RSA_Data.PRT_ROI = tmp_strct.(['RSA_Data_',Subj_names{sub}]).PRT.curROI;
    else
        RSA_Data.Encoding_Data{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.Encoding_Data;
        RSA_Data.TrialInfo{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;
    end
    
    if(~isempty(tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full))
        RSA_Data.OCC.full_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full{1,msr},[3 1 2]);
        RSA_Data.TMP.full_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_full{1,msr},[3 1 2]);
        RSA_Data.FRT.full_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_full{1,msr},[3 1 2]);
        RSA_Data.CNT.full_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).CNT.RSA_full{1,msr},[3 1 2]);
        RSA_Data.PRT.full_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_full{1,msr},[3 1 2]);
    end
    RSA_Data.OCC.red16_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_16{1,msr},[3 1 2]);
    RSA_Data.TMP.red16_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_16{1,msr},[3 1 2]);
    RSA_Data.FRT.red16_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_16{1,msr},[3 1 2]);
    RSA_Data.CNT.red16_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).CNT.RSA_16{1,msr},[3 1 2]);
    RSA_Data.PRT.red16_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_16{1,msr},[3 1 2]);
end

    
    
%% Hypotheses Matrix
trl_mat = [kron([1;2],ones(64,1)) kron([1;2;1;2],ones(32,1))];
Perceptual_Mat_full = zeros(size(trl_mat,1));
Semantic_Mat_full = zeros(size(trl_mat,1));
for i = 1:size(trl_mat,1)-1
    for j = (i+1):size(trl_mat,1)
        if(j ~= size(trl_mat,1) - (i - 1))
            if(trl_mat(i,1) == 1 && trl_mat(j,1) == 1)
                Perceptual_Mat_full(i,j) = 1;
            elseif(trl_mat(i,1) == 2 && trl_mat(j,1) == 2)
                Perceptual_Mat_full(i,j) = 2;
            else
                Perceptual_Mat_full(i,j) = -1;
            end

            if(trl_mat(i,2) == 1 && trl_mat(j,2) == 1)
                Semantic_Mat_full(i,j) = 1;
            elseif(trl_mat(i,2) == 2 && trl_mat(j,2) == 2)
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

Perceptual_Mat_red16 = zeros(16);
Semantic_Mat_red16 = zeros(16);
for i = 1:16-1
    for j = (i+1):16
        if(j ~= 16 - (i - 1))
            if(trl_mat(8*i,1) == 1 && trl_mat(8*j,1) == 1)
                Perceptual_Mat_red16(i,j) = 1;
            elseif(trl_mat(8*i,1) == 2 && trl_mat(8*j,1) == 2)
                Perceptual_Mat_red16(i,j) = 2;
            else
                Perceptual_Mat_red16(i,j) = -1;
            end

            if(trl_mat(8*i,2) == 1 && trl_mat(8*j,2) == 1)
                Semantic_Mat_red16(i,j) = 1;
            elseif(trl_mat(8*i,2) == 2 && trl_mat(8*j,2) == 2)
                Semantic_Mat_red16(i,j) = 2;
            else
                Semantic_Mat_red16(i,j) = -1;
            end
        end
    end
end

figure
subplot(1,2,1)
imagesc(Perceptual_Mat_red16); title('Perceptual Hypothesis Matrix')
axis square
subplot(1,2,2)
imagesc(Semantic_Mat_red16); title('Semantic Hypothesis Matrix')
axis square

 
    
%% Create RSA Time Courses

RSA_Time = [];
TimeVec = RSA_Data.TimeVec;
for sub = 1:length(Subj_names)
    
    cur_data = zeros(size(RSA_Data.OCC.red16_Data,3));
    for tp = 1:length(TimeVec)
        
        % Occipital
        cur_data = squeeze(RSA_Data.OCC.red16_Data(sub,tp,:,:));
        % Perceptual Dimension
        RSA_Time.OCC.Perceptual.Drawing(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 == 1));
        RSA_Time.OCC.Perceptual.Picture(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 == 2));
        RSA_Time.OCC.Perceptual.Within(sub,tp)      = nanmean(cur_data(Perceptual_Mat_red16 > 0));
        RSA_Time.OCC.Perceptual.Between(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 < 0));
        % Semantic Dimension
        RSA_Time.OCC.Semantic.Animate(sub,tp)       = nanmean(cur_data(Semantic_Mat_red16 == 1));
        RSA_Time.OCC.Semantic.Inanimate(sub,tp)     = nanmean(cur_data(Semantic_Mat_red16 == 2));
        RSA_Time.OCC.Semantic.Within(sub,tp)        = nanmean(cur_data(Semantic_Mat_red16 > 0));
        RSA_Time.OCC.Semantic.Between(sub,tp)       = nanmean(cur_data(Semantic_Mat_red16 < 0));
        
        % Temporal
        cur_data = squeeze(RSA_Data.TMP.red16_Data(sub,tp,:,:));
        % Perceptual Dimension
        RSA_Time.TMP.Perceptual.Drawing(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 == 1));
        RSA_Time.TMP.Perceptual.Picture(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 == 2));
        RSA_Time.TMP.Perceptual.Within(sub,tp)      = nanmean(cur_data(Perceptual_Mat_red16 > 0));
        RSA_Time.TMP.Perceptual.Between(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 < 0));
        % Semantic Dimension
        RSA_Time.TMP.Semantic.Animate(sub,tp)       = nanmean(cur_data(Semantic_Mat_red16 == 1));
        RSA_Time.TMP.Semantic.Inanimate(sub,tp)     = nanmean(cur_data(Semantic_Mat_red16 == 2));
        RSA_Time.TMP.Semantic.Within(sub,tp)        = nanmean(cur_data(Semantic_Mat_red16 > 0));
        RSA_Time.TMP.Semantic.Between(sub,tp)       = nanmean(cur_data(Semantic_Mat_red16 < 0));
        
        % Frontal
        cur_data = squeeze(RSA_Data.FRT.red16_Data(sub,tp,:,:));
        % Perceptual Dimension
        RSA_Time.FRT.Perceptual.Drawing(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 == 1));
        RSA_Time.FRT.Perceptual.Picture(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 == 2));
        RSA_Time.FRT.Perceptual.Within(sub,tp)      = nanmean(cur_data(Perceptual_Mat_red16 > 0));
        RSA_Time.FRT.Perceptual.Between(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 < 0));
        % Semantic Dimension
        RSA_Time.FRT.Semantic.Animate(sub,tp)       = nanmean(cur_data(Semantic_Mat_red16 == 1));
        RSA_Time.FRT.Semantic.Inanimate(sub,tp)     = nanmean(cur_data(Semantic_Mat_red16 == 2));
        RSA_Time.FRT.Semantic.Within(sub,tp)        = nanmean(cur_data(Semantic_Mat_red16 > 0));
        RSA_Time.FRT.Semantic.Between(sub,tp)       = nanmean(cur_data(Semantic_Mat_red16 < 0));
        
        % Central
        cur_data = squeeze(RSA_Data.CNT.red16_Data(sub,tp,:,:));
        % Perceptual Dimension
        RSA_Time.CNT.Perceptual.Drawing(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 == 1));
        RSA_Time.CNT.Perceptual.Picture(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 == 2));
        RSA_Time.CNT.Perceptual.Within(sub,tp)      = nanmean(cur_data(Perceptual_Mat_red16 > 0));
        RSA_Time.CNT.Perceptual.Between(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 < 0));
        % Semantic Dimension
        RSA_Time.CNT.Semantic.Animate(sub,tp)       = nanmean(cur_data(Semantic_Mat_red16 == 1));
        RSA_Time.CNT.Semantic.Inanimate(sub,tp)     = nanmean(cur_data(Semantic_Mat_red16 == 2));
        RSA_Time.CNT.Semantic.Within(sub,tp)        = nanmean(cur_data(Semantic_Mat_red16 > 0));
        RSA_Time.CNT.Semantic.Between(sub,tp)       = nanmean(cur_data(Semantic_Mat_red16 < 0));
        
        % Parietal
        cur_data = squeeze(RSA_Data.PRT.red16_Data(sub,tp,:,:));
        % Perceptual Dimension
        RSA_Time.PRT.Perceptual.Drawing(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 == 1));
        RSA_Time.PRT.Perceptual.Picture(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 == 2));
        RSA_Time.PRT.Perceptual.Within(sub,tp)      = nanmean(cur_data(Perceptual_Mat_red16 > 0));
        RSA_Time.PRT.Perceptual.Between(sub,tp)     = nanmean(cur_data(Perceptual_Mat_red16 < 0));
        % Semantic Dimension
        RSA_Time.PRT.Semantic.Animate(sub,tp)       = nanmean(cur_data(Semantic_Mat_red16 == 1));
        RSA_Time.PRT.Semantic.Inanimate(sub,tp)     = nanmean(cur_data(Semantic_Mat_red16 == 2));
        RSA_Time.PRT.Semantic.Within(sub,tp)        = nanmean(cur_data(Semantic_Mat_red16 > 0));
        RSA_Time.PRT.Semantic.Between(sub,tp)       = nanmean(cur_data(Semantic_Mat_red16 < 0));
        
    end

end


%% Permutation Matrix


cfg = [];
cfg.nPerms = 1000;
cfg.thresh_pval = 0.01;
cfg.mcc_cluster_pval = 0.01;
cfg.TimeVec = RSA_Data.TimeVec;
cfg.Hyp_Mat = double(Semantic_Mat_red16>0) + -double(Semantic_Mat_red16<0);
cfg.matshuff = false;
cfg.twoside = true;
%Results = rsa_perm(cfg, RSA_Data.TMP.red16_Data);


% Plot RSA Time series

ROI = {'OCC','TMP','FRT','CNT','PRT'}; r = 3;
Cat = {'Perceptual','Semantic'}; c = 1;
dt = 3:4;

figure('Pos', [325 510 650 402]);
plot([]); hold on
Dat_names = fieldnames(RSA_Time.(ROI{r}).(Cat{c}));
dat1 = nanmean(RSA_Time.(ROI{r}).(Cat{c}).(Dat_names{dt(1)})(:,:),1);
dat2 = nanmean(RSA_Time.(ROI{r}).(Cat{c}).(Dat_names{dt(2)})(:,:),1);
SEM1 = nanstd(RSA_Time.(ROI{r}).(Cat{c}).(Dat_names{dt(1)})(:,:),0,1)./sqrt(size(RSA_Time.(ROI{r}).(Cat{c}).(Dat_names{dt(1)}),1));
SEM2 = nanstd(RSA_Time.(ROI{r}).(Cat{c}).(Dat_names{dt(1)})(:,:),0,1)./sqrt(size(RSA_Time.(ROI{r}).(Cat{c}).(Dat_names{dt(1)}),1));
fill([TimeVec fliplr(TimeVec)],[dat1 fliplr(dat1 + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[dat1 fliplr(dat1 - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[dat2 fliplr(dat2 + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[dat2 fliplr(dat2 - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
h1 = plot(TimeVec, dat1,'b','linewidth',2);
h2 = plot(TimeVec, dat2,'r','linewidth',2);
hold off
ylabel(RSA_Data.meas16{msr}); xlabel('Time (s)'); title([Cat{c},': ',Dat_names{dt(1)},' vs ',Dat_names{dt(2)}])
xlim([-0.2 1.5]);lg = legend([h1 h2], {Dat_names{dt(1)},Dat_names{dt(2)}}); legend boxoff; set(lg,'FontSize',14)
box off; %ylim([0.42 0.6])
set(gca,'linewidth',2.5,'FontSize',14)
%saveas(gcf,sprintf('Results/%s_16_%s_LDA_%s.png',ROI{r},Cat{c},datanames{d}))
%close(gcf)


% Plot RSA Difference Time series

r1 = 3; r2 = 5;
c1 = 1; c2 = 2;
dt1 = [3 4]; dt2 = [3 4];

if(sum(dt1) == 7)
    Hyp_Mat_per = double(Perceptual_Mat_red16>0) + -double(Perceptual_Mat_red16<0);
    Hyp_Mat_sem = double(Semantic_Mat_red16>0) + -double(Semantic_Mat_red16<0);
else
    Hyp_Mat_per = double(Perceptual_Mat_red16==1) + -double(Perceptual_Mat_red16==2);
    Hyp_Mat_sem = double(Semantic_Mat_red16==1) + -double(Semantic_Mat_red16==2);
end

cfg = [];
cfg.nPerms = 1000;
cfg.thresh_pval = 0.05;
cfg.mcc_cluster_pval = 0.05;
cfg.TimeVec = RSA_Data.TimeVec;
cfg.Hyp_Mat = Hyp_Mat_per;
cfg.matshuff = false;
cfg.twoside = true;
Results1 = rsa_perm(cfg, RSA_Data.(ROI{r1}).red16_Data);

cfg.Hyp_Mat = Hyp_Mat_sem;
Results2 = rsa_perm(cfg, RSA_Data.(ROI{r2}).red16_Data);

figure('Pos', [325 510 650 402]);
plot([]); hold on
Dat_names1 = fieldnames(RSA_Time.(ROI{r1}).(Cat{c1}));
dat1 = RSA_Time.(ROI{r1}).(Cat{c1}).(Dat_names1{dt1(2)}) - RSA_Time.(ROI{r1}).(Cat{c1}).(Dat_names1{dt1(1)});
Dat_names2 = fieldnames(RSA_Time.(ROI{r2}).(Cat{c2}));
dat2 = RSA_Time.(ROI{r2}).(Cat{c2}).(Dat_names2{dt2(2)}) - RSA_Time.(ROI{r2}).(Cat{c2}).(Dat_names2{dt2(1)});
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
%saveas(gcf,'Results/Enc_TMP_PerceptualvsSemantic_BT-WI_Dat16_LDA.png')
%close(gcf)



%% Plot RSA Matrices Timeseries


ROI = {'OCC','TMP','FRT','CNT','PRT'}; 

[a,b] = unique(RSA_Data.TrialInfo{1,1}(1:64,10),'stable');
Cat_names = strcat(strrep(strrep(RSA_Data.TrialInfo{1,1}([b;b+64],6),'Picture','Pic'),'Drawing','Dr'), {' - '}, [a;a]);


% Plot Single GA Matrix

r = 1;
tp = 0.2; % s
clim_wi = true;

if(~clim_wi)
    calc_dat = nanmean(RSA_Data.(ROI{r}).red16_Data(:,:),1);
    colb_lim = prctile(calc_dat(calc_dat ~= 0), [2.5 97.5]);
end

figure
[~,t_idx] = min(abs(RSA_Data.TimeVec - tp));
cur_data = squeeze(nanmean(RSA_Data.(ROI{r}).red16_Data(:,t_idx,:,:),1));
if(clim_wi)
    colb_lim = prctile(cur_data(cur_data(:) ~= 0), [2.5 97.5]);
end
imagesc(cur_data+cur_data');
colorbar
axis square
set(gca,'clim',colb_lim,'xtick',1:16,'xticklabels',Cat_names,'XTickLabelRotation',55,...
    'ytick',1:16,'yticklabels',Cat_names, 'TickLength',[0 0])



% Plot Timeseries GA Matrix

r = 2;

if(~clim_wi)
    calc_dat = nanmean(RSA_Data.(ROI{r}).red16_Data(:,:),1);
    colb_lim = prctile(calc_dat(calc_dat ~= 0), [2.5 97.5]);
end

figure('pos',[202 73 1481 898])
tps = -0.2:0.1:1.2; %s
for sbp = 1:length(tps)
    [~,t_idx] = min(abs(RSA_Data.TimeVec - tps(sbp)));
    cur_data = squeeze(nanmean(RSA_Data.(ROI{r}).red16_Data(:,t_idx,:,:),1));
    if(clim_wi)
        colb_lim = prctile(cur_data(cur_data(:) ~= 0), [2.5 97.5]);
    end
    subplot(3,5,sbp);
    imagesc(cur_data+cur_data');
    axis square
    set(gca,'clim',colb_lim,'xtick',[],'ytick',[],'TickLength',[0 0])
    if(sbp == 11)
        set(gca,'xtick',1:16,'xticklabels',Cat_names,'XTickLabelRotation',55,...
            'ytick',1:16,'yticklabels',Cat_names, 'TickLength',[0 0])
    end
    title(sprintf('%4.0f ms', tps(sbp)*1000))
end
cl = colorbar; set(cl,'Pos', [0.9 0.2 0.1 0.1])
pos_val = get(gca,'Position');
set(cl,'Pos', [0.92 pos_val(2) 0.0119 pos_val(4)],'Xtick',linspace(colb_lim(1),colb_lim(2),5))




%% Plot Shuffled Hypotheses Matrices
% Test_Mat = Semantic_Mat_16;
% %Test_Mat = logical(Perceptual_Mat_16).*(rand(16)*1000);
% %labels = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p'};
% %Test_Mat = [0 1 2 3 4 0; 0 0 5 6 0 7; 0 0 0 0 8 9; 0 0 0 0 10 11; 0 0 0 0 0 12; 0 0 0 0 0 0];
% %rand_idx = [1 4 3 5 2 6];
% 
% figure
% subplot(4,4,1)
% imagesc(Test_Mat)
% axis square
% for sbp = 2:16
% 
%     Shuff_Mat = zeros(size(Test_Mat));
%     rand_idx = randperm(size(Test_Mat,1));
%     for row = 1:size(Test_Mat,1)-1
%         for col = (row+1):size(Test_Mat,1)
%             if(Test_Mat(rand_idx(row),rand_idx(col)) ~= 0)
%                 Shuff_Mat(row,col) = Test_Mat(rand_idx(row),rand_idx(col));
%             else
%                 Shuff_Mat(row,col) = Test_Mat(rand_idx(col),rand_idx(row));
%             end
%         end
%     end
%     
%    
%     subplot(4,4,sbp)
%     imagesc(Shuff_Mat)
%     axis square
% end
% 
% lab_idx = [find(strcmp(labels,'d')), find(strcmp(labels,'a'))];
% Test_Mat(min(lab_idx), max(lab_idx))
% lab_perm_idx = [find(strcmp(labels(rand_idx),'d')), find(strcmp(labels(rand_idx),'a'))];
% Shuff_Mat(min(lab_perm_idx), max(lab_perm_idx))