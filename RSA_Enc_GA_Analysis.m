%% RSA Encoding - All Subjects Analysis

%% Subject Names
load('RSA_Data_Enc','Subj_names')


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs
ROI_idx = {ROI_occ_idx; ROI_temp_idx; ROI_front_idx; ROI_pari_idx};
%ROI_idx = {ROI_all_idx};

%% Hypothesis Matrix
Hypothesis_Matrix


%% Measure
msr = 1;


%% Matrix Indices to choose
UTmat = triu(ones(16)).*~eye(16);
UTM_idx = find(UTmat(:) == 1);
SensUTmat = UTmat(1:end-1,:);
SensUTM_idx = find(SensUTmat(:) == 1);

%% Create Data
RSA_Data_Enc = [];
for sub = 1:length(Subj_names)
    tmp_strct = load(['RSA_Data/RSA_Data_Enc_',Subj_names{sub}]);
    if(sub == 1)
        ROI = fieldnames(tmp_strct.(['RSA_Data_',Subj_names{1}]));
        RSA_Data_Enc.Names   = Subj_names;
        for fn = fieldnames(tmp_strct.(['RSA_Data_',Subj_names{sub}]).(ROI{1}))'
            RSA_Data_Enc.(fn{1}) = tmp_strct.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).(fn{1});
        end
        RSA_Data_Enc = rmfield(RSA_Data_Enc, {'Name','TimeVec1024','RSA_full','MDS_full','RSA_16','Encoding_Data','TrialInfo'...
                                      'rsa_dim','curROI','curROI_name'});
        RSA_Data_Enc.Encoding_Data{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).Encoding_Data;
        RSA_Data_Enc.TrialInfo{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).TrialInfo;
        [~,RSA_Data_Enc.Cond_idx.Perc_WI_BT_idx(:,1),~] = intersect(UTM_idx,find(Perceptual_Mat_red16(:) > 0));
        [~,RSA_Data_Enc.Cond_idx.Perc_WI_BT_idx(:,2),~] = intersect(UTM_idx,find(Perceptual_Mat_red16(:) < 0));
        [~,RSA_Data_Enc.Cond_idx.Sem_WI_BT_idx(:,1),~]  = intersect(UTM_idx,find(Semantic_Mat_red16(:) > 0));
        [~,RSA_Data_Enc.Cond_idx.Sem_WI_BT_idx(:,2),~]  = intersect(UTM_idx,find(Semantic_Mat_red16(:) < 0));
        [~,RSA_Data_Enc.Cond_idx.Perc_Dr_Ph_idx(:,1),~] = intersect(UTM_idx,find(Perceptual_Mat_red16(:) == 1));
        [~,RSA_Data_Enc.Cond_idx.Perc_Dr_Ph_idx(:,2),~] = intersect(UTM_idx,find(Perceptual_Mat_red16(:) == 2));
        [~,RSA_Data_Enc.Cond_idx.Sem_An_Ia_idx(:,1),~]  = intersect(UTM_idx,find(Semantic_Mat_red16(:) == 1));
        [~,RSA_Data_Enc.Cond_idx.Sem_An_Ia_idx(:,2),~]  = intersect(UTM_idx,find(Semantic_Mat_red16(:) == 2));
        for r = 1:length(ROI)
            RSA_Data_Enc.([ROI{r},'_ROI']) = tmp_strct.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).curROI;
        end
    else
        RSA_Data_Enc.Encoding_Data{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).Encoding_Data;
        RSA_Data_Enc.TrialInfo{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).TrialInfo;
    end
    
    for r = 1:length(ROI)
        if(~isempty(tmp_strct.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_full))
            hlp = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_full{1,msr},[3 1 2]);
            RSA_Data_Enc.(ROI{r}).full_Data(sub,:,:) = hlp(:,UTM_idx);
        end
        hlp = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_16{1,msr},[3 1 2]);
        RSA_Data_Enc.(ROI{r}).red16_Data(sub,:,:) = hlp(:,UTM_idx);
        if(msr == 1)
            hlp = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_16{3,msr},[4 3 1 2]);
            hlp2 = hlp(:,:,SensUTM_idx); %bsxfun(@times, hlp(:,:,SensUTM_idx),RSA_Data.(ROI{r}).red16_Data(sub,:,:));
            RSA_Data_Enc.(ROI{r}).red16_SensorData(sub,:,:,1) = mean(hlp2(:,:,RSA_Data_Enc.Cond_idx.Perc_WI_BT_idx(:,1)),3); %bsxfun(@times, mean(hlp2(:,:,RSA_Data.Cond_idx.Perc_WI_BT_idx(:,1)),3), mean(RSA_Data.(ROI{r}).red16_Data(sub,:,RSA_Data.Cond_idx.Perc_WI_BT_idx(:,2)),3) - mean(RSA_Data.(ROI{r}).red16_Data(sub,:,RSA_Data.Cond_idx.Perc_WI_BT_idx(:,1)),3));
            RSA_Data_Enc.(ROI{r}).red16_SensorData(sub,:,:,2) = mean(hlp2(:,:,RSA_Data_Enc.Cond_idx.Perc_WI_BT_idx(:,2)),3); %bsxfun(@times, mean(hlp2(:,:,RSA_Data.Cond_idx.Perc_WI_BT_idx(:,2)),3), mean(RSA_Data.(ROI{r}).red16_Data(sub,:,RSA_Data.Cond_idx.Perc_WI_BT_idx(:,2)),3) - mean(RSA_Data.(ROI{r}).red16_Data(sub,:,RSA_Data.Cond_idx.Perc_WI_BT_idx(:,1)),3));
            RSA_Data_Enc.(ROI{r}).red16_SensorData(sub,:,:,3) = mean(hlp2(:,:,RSA_Data_Enc.Cond_idx.Sem_WI_BT_idx(:,1)),3); %bsxfun(@times, mean(hlp2(:,:,RSA_Data.Cond_idx.Sem_WI_BT_idx(:,1)),3), mean(RSA_Data.(ROI{r}).red16_Data(sub,:,RSA_Data.Cond_idx.Sem_WI_BT_idx(:,2)),3) - mean(RSA_Data.(ROI{r}).red16_Data(sub,:,RSA_Data.Cond_idx.Sem_WI_BT_idx(:,1)),3));
            RSA_Data_Enc.(ROI{r}).red16_SensorData(sub,:,:,4) = mean(hlp2(:,:,RSA_Data_Enc.Cond_idx.Sem_WI_BT_idx(:,2)),3); %bsxfun(@times, mean(hlp2(:,:,RSA_Data.Cond_idx.Sem_WI_BT_idx(:,2)),3), mean(RSA_Data.(ROI{r}).red16_Data(sub,:,RSA_Data.Cond_idx.Sem_WI_BT_idx(:,2)),3) - mean(RSA_Data.(ROI{r}).red16_Data(sub,:,RSA_Data.Cond_idx.Sem_WI_BT_idx(:,1)),3));
        end
    end
    clear tmp_strct
end

    



%% Plotting Results

% Parameters
TimeVec_Enc = RSA_Data_Enc.TimeVec;
Cat = {'Perc','Sem'};
Comp = {'WI_BT','Dr_Ph'; 'WI_BT','An_Ia'};
%ROI_names = {'Occipital','Temporal','Frontal','Parietal'};
ROI_names = {'Wholehead'};
Comp_names = {{'Within','Between'}, {'Drawing','Photograph'}; {'Within','Between'}, {'Animate','Inanimate'}};
Cat_names = {'Perceptual','Semantic'};



%% Plot RSA Time series
r = 4;
c = 1;
cp = 1;

figure('Pos', [325 210 650 402]);
plot([]); hold on
dat1 = mean(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,1)),3),1);
dat2 = mean(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,2)),3),1);
SEM1 = nanstd(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,1)),3),0,1)./sqrt(size(RSA_Data_Enc.(ROI{r}).red16_Data,1));
SEM2 = nanstd(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,2)),3),0,1)./sqrt(size(RSA_Data_Enc.(ROI{r}).red16_Data,1));
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat1 fliplr(dat1 + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat1 fliplr(dat1 - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat2 fliplr(dat2 + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat2 fliplr(dat2 - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
h1 = plot(TimeVec_Enc, dat1,'b','linewidth',2);
h2 = plot(TimeVec_Enc, dat2,'r','linewidth',2);
hold off
ylabel(RSA_Data_Enc.meas16{msr}); xlabel('Time (s)'); %title([ROI{r}])
xlim([-0.2 1.5]);lg = legend([h1 h2], {[Cat_names{c},' ',Comp_names{c,cp}{1}],[Cat_names{c},' ',Comp_names{c,cp}{2}]}); legend boxoff; set(lg,'FontSize',14)
box off; %ylim([0.42 0.6])
set(gca,'linewidth',2.5,'FontSize',14)
%saveas(gcf,sprintf('Results/Enc_%s_16_%s_LDA_%s.png',ROI{r},Cat{c},datanames{d}))
%close(gcf)



%% Plot RSA Difference Time series
r1 = 4; r2 = 4;
c1 = 1; c2 = 2;
cp = 1;

if(cp == 1)
    Perceptual_idx  = RSA_Data_Enc.Cond_idx.Perc_WI_BT_idx;
    Semantic_idx    = RSA_Data_Enc.Cond_idx.Sem_WI_BT_idx;
else
    Perceptual_idx  = RSA_Data_Enc.Cond_idx.Perc_Dr_Ph_idx;
    Semantic_idx 	= RSA_Data_Enc.Cond_idx.Sem_An_Ia_idx;
end

cfg = [];
cfg.nPerms = 1000;
cfg.thresh_pval = 0.05;
cfg.mcc_cluster_pval = 0.05;
cfg.TimeVec = RSA_Data_Enc.TimeVec;
cfg.Hyp_idx = Perceptual_idx;
cfg.matshuff = false;
cfg.twoside = true;
Results1 = rsa_perm(cfg, RSA_Data_Enc.(ROI{r1}).red16_Data);

cfg.Hyp_idx = Semantic_idx;
Results2 = rsa_perm(cfg, RSA_Data_Enc.(ROI{r2}).red16_Data);

figure('Pos', [114 259 1261 714]); %[325 510 650 402]
plot([]); hold on
dat1 = mean(RSA_Data_Enc.(ROI{r1}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c1},'_',Comp{cp},'_idx'])(:,2)),3) - mean(RSA_Data_Enc.(ROI{r1}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c1},'_',Comp{cp},'_idx'])(:,1)),3);
dat2 = mean(RSA_Data_Enc.(ROI{r2}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c2},'_',Comp{cp},'_idx'])(:,2)),3) - mean(RSA_Data_Enc.(ROI{r2}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c2},'_',Comp{cp},'_idx'])(:,1)),3);
SEM1 = nanstd(dat1,0,1)./sqrt(size(dat1,1)); SEM2 = nanstd(dat2,0,1)./sqrt(size(dat2,1));
% SEM1 = zeros(1,length(TimeVec)); SEM2 = zeros(1,length(TimeVec));
% for tp = 1:length(TimeVec)
%     [~,~,~,stats] = ttest(RSA_Time.(ROI{r1}).(Cat{c1}).(Dat_names{dt(2)})(:,tp), RSA_Time.(ROI{r1}).(Cat{c1}).(Dat_names{dt(1)})(:,tp));
%     SEM1(1,tp) = tinv(.95,stats.df)*(stats.sd/sqrt(size(dat1,1)));
%     [~,~,~,stats] = ttest(RSA_Time.(ROI{r2}).(Cat{c2}).(Dat_names{dt(2)})(:,tp), RSA_Time.(ROI{r2}).(Cat{c2}).(Dat_names{dt(1)})(:,tp));
%     SEM2(1,tp) = tinv(.95,stats.df)*(stats.sd/sqrt(size(dat2,1)));
% end
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
xlim([-0.2 1.5]); %ylim([0.35 0.6]);
plot([-0.2 1.5],[0 0],'--k','linewidth',1)
h1 = plot(TimeVec_Enc, nanmean(dat1,1),'b','linewidth',2);
h2 = plot(TimeVec_Enc, nanmean(dat2,1),'r','linewidth',2);
ylabel(['Diff ',RSA_Data_Enc.meas16{msr}]); xlabel('Time (s)'); 
lg = legend([h1 h2], {[ROI{r1},' ',Cat{c1},' ',strjoin(fliplr(Comp_names{c,cp}),' - ')],[ROI{r2},' ',Cat{c2},' ',strjoin(fliplr(Comp_names{c,cp}),' - ')]}); legend boxoff; set(lg,'FontSize',14)
box off;
set(gca,'linewidth',2.5,'FontSize',14,'xlim',[-0.2 1],'xtick',[-0.2:0.1:1])
if(Results1.H == 1)
    sign_mcc_clust = abs(Results1.zmapthresh);
    sign_mcc_clust(sign_mcc_clust > 0) = min(get(gca,'ylim'))*0.8;
    plot(TimeVec_Enc,sign_mcc_clust,'bo','MarkerFaceColor','b')
end
if(Results2.H == 1)
    sign_mcc_clust = abs(Results2.zmapthresh);
    sign_mcc_clust(sign_mcc_clust > 0) = min(get(gca,'ylim'))*0.9;
    plot(TimeVec_Enc,sign_mcc_clust,'ro','MarkerFaceColor','r')
end
hold off
%saveas(gcf,'Results/Enc_OCC_TMP_PerceptualvsSemantic_BT-WI_Dat16_eucl.png')
%close(gcf)



%% Plot RSA Matrices Timeseries


ROI = {'OCC','TMP','FRT','PRT'}; 

[a,b] = unique(RSA_Data_Enc.TrialInfo{1,1}(1:64,10),'stable');
Cat_names = strcat(strrep(strrep(RSA_Data_Enc.TrialInfo{1,1}([b;b+64],6),'Picture','Pic'),'Drawing','Dr'), {' - '}, [a;a]);


% Plot Single GA Matrix

r = 4;
tp = 0.2; % s

calc_dat = nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:),1);
calc_dat(calc_dat == 0) = [];
colb_lim = prctile(calc_dat, [2.5 97.5]);

figure
cur_data = zeros(16);
cur_data(UTM_idx) = squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,dsearchn(RSA_Data_Enc.TimeVec', tp),:),1));
imagesc(cur_data);
colorbar
axis square
set(gca,'clim',colb_lim,'xtick',1:16,'xticklabels',Cat_names,'XTickLabelRotation',55,...
    'ytick',1:16,'yticklabels',Cat_names, 'TickLength',[0 0])
title('OCC Encoding at 200 ms')






ROI = {'OCC','TMP','FRT','PRT'}; 

[a,b] = unique(RSA_Data_Enc.TrialInfo{1,1}(1:64,10),'stable');
Cat_names = strcat(strrep(strrep(RSA_Data_Enc.TrialInfo{1,1}([b;b+64],6),'Picture','Pic'),'Drawing','Dr'), {' - '}, [a;a]);


% Plot Single GA Matrix

r = 1;
tp = 0.2; % s
clim_wi = true;

if(~clim_wi)
    calc_dat = nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:),1);
    colb_lim = prctile(calc_dat(calc_dat ~= 0), [2.5 97.5]);
end

figure
cur_data = zeros(16);
cur_data(UTM_idx) = squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,dsearchn(RSA_Data_Enc.TimeVec', tp),:),1));
if(clim_wi)
    colb_lim = prctile(cur_data(cur_data(:) ~= 0), [2.5 97.5]);
end
imagesc(cur_data+cur_data');
colorbar
axis square
set(gca,'clim',colb_lim,'xtick',1:16,'xticklabels',Cat_names,'XTickLabelRotation',55,...
    'ytick',1:16,'yticklabels',Cat_names, 'TickLength',[0 0])
title(sprintf('%s Encoding at %4.0f ms',ROI{r},round(1000*tp)))


% Plot Timeseries GA Matrix

r = 2;

if(~clim_wi)
    calc_dat = nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:),1);
    colb_lim = prctile(calc_dat(calc_dat ~= 0), [2.5 97.5]);
end

figure('pos',[202 73 1481 898])
tps = -0.2:0.1:1.2; %s
for sbp = 1:length(tps)
    cur_data = zeros(16);
    cur_data(UTM_idx) = squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,dsearchn(RSA_Data_Enc.TimeVec', tps(sbp)),:,:),1));
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


% Simulate RDMs for SFN Poster
for i = 1:6
    figure('pos',[200 200 400 400])
    cur_data = (0.4*rand(16)+0.4).*(~eye(16));
    cur_data = triu(cur_data) + triu(cur_data)';
    imagesc(cur_data)
    set(gca,'Pos', [0 0 1 1])
    axis square
    set(gca,'clim',[0.3 0.9],'xtick','','ytick','','TickLength',[0 0])
    saveas(gcf,sprintf('Results/SFN_Poster/Ret_RDM_plain_Ex_%i.png',i))
    close(gcf)
end

figure('pos',[834   697   606    32])
cur_data = (0.4*rand(16)+0.4).*(~eye(16));
cur_data = triu(cur_data) + triu(cur_data)';
imagesc(cur_data(:,end)')
set(gca,'xtick',[],'ytick',[])

cx = 0; cy = 0;
for i = 1:2
    figure('pos',[200 200 400 400])
    axes()
    set(gca,'Pos', [0 0 1 1])
    axis square
    set(gca,'clim',[0.1 0.9],'xtick','','ytick','','TickLength',[0 0])
    hold on
    if(i == 1)
        rectangle('Position',[8.55 0.55 7.9 7.9],'EdgeColor',[1 0 0],'LineWidth',4)
        fill([8.55 8.55+7.9 8.55+7.9 8.55],[0.55 0.55 0.55+7.9 0.55+7.9],'red','FaceAlpha',0.3)
        rectangle('Position',[0.55 8.55 7.9 7.9],'EdgeColor',[1 0 0],'LineWidth',4)
        fill([0.55 0.55+7.9 0.55+7.9 0.55],[8.55 8.55 8.55+7.9 8.55+7.9],'red','FaceAlpha',0.3)
        rectangle('Position',[0.55 0.55 7.9 7.9],'EdgeColor',[0 0 1],'LineWidth',4)
        fill([0.55 0.55+7.9 0.55+7.9 0.55],[0.55 0.55 0.55+7.9 0.55+7.9],'blue','FaceAlpha',0.3)
        rectangle('Position',[8.55 8.55 7.9 7.9],'EdgeColor',[0 0 1],'LineWidth',4)
        fill([8.55 8.55+7.9 8.55+7.9 8.55],[8.55 8.55 8.55+7.9 8.55+7.9],'blue','FaceAlpha',0.3)
        hold off
        saveas(gcf,sprintf('Results/SFN_Poster/Ret_RDM_btwi_Perc_Ex_%i.png',i))
        close(gcf)
    else
        for j = 1:16
            px = 4*cx + 0.55; py = 4*cy + 0.55;
            rectangle('Position',[ px  py 3.9 3.9],'EdgeColor',[mod(cx+cy,2) 0 ~mod(cx+cy,2)],'LineWidth',4)
            fill([px px+3.9 px+3.9 px],[py py py+3.9 py+3.9],[mod(cx+cy,2) 0 ~mod(cx+cy,2)],'FaceAlpha',0.3)
            cx = mod(cx + 1,4);
            cy = cy + (cx == 0);
        end
        hold off
        saveas(gcf,sprintf('Results/SFN_Poster/Ret_RDM_btwi_Sem_Ex_%i.png',i-1))
        close(gcf)
    end
end




%% RDM with Categories


ROI = {'OCC','TMP','FRT','PRT'};
Dim_names = {};

[a,b] = unique(RSA_Data_Enc.TrialInfo{1,1}(1:64,6),'stable');
Dim_names{1,1} = RSA_Data_Enc.TrialInfo{1,1}([b;b+64],6);
Dim_names{1,2} = [8.5; 8.5];
[a,b] = unique(RSA_Data_Enc.TrialInfo{1,1}(1:64,8),'stable');
Dim_names{2,1} = RSA_Data_Enc.TrialInfo{1,1}([b;b+64],8);
Dim_names{2,2} = [4.5:4:12.5; 4.5:4:12.5];


% Plot Single GA Matrix
r = 1;
c = 1;
tp = 0.2; % s
clim_wi = true;
if(~clim_wi)
    calc_dat = nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:),1);
    colb_lim = prctile(calc_dat(calc_dat ~= 0), [2.5 97.5]);
end

figure('Units','pixels')
h1 = axes('Units','pixels','xlim',[0 18],'ylim',[0 18]);
cur_data = zeros(16);
cur_data(UTM_idx) = squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,dsearchn(RSA_Data_Enc.TimeVec', tp),:),1));
if(clim_wi)
    colb_lim = [min(cur_data(cur_data(:) ~= 0))*.95 prctile(cur_data(cur_data(:) ~= 0), 90)];
end
imagesc(cur_data);
axis square
grid on
set(get(gca,'Yruler'),'Minortick','off') %Dim_names{2,2}(1,:))
set(get(gca,'Xruler'),'Minortick','off')%Dim_names{2,2}(1,:))
set(gca,'clim',colb_lim,'xtick',Dim_names{1,2}(1,:),'xticklabels',[],...
        'ytick',Dim_names{1,2}(2,:),'yticklabels',[], 'TickLength',[0 0],'XMinorgrid','on','YMinorgrid','on',...
        'gridcolor','w','gridalpha',.9,'minorgridlinestyle','--','minorgridalpha',.5,'MinorGridColor','w')
pos = get(h1,'pos') + [-40 -40 0 0];
set(h1,'pos',pos)
text(3,-1.4,Dim_names{1,1}{1},'fontsize',13)
text(11,-1.4,Dim_names{1,1}{2},'fontsize',13)
text(18.2,3,Dim_names{1,1}{1},'fontsize',13,'Rotation',-90)
text(18.2,11,Dim_names{1,1}{2},'fontsize',13,'Rotation',-90)
text(1,-0.3,Dim_names{2,1}{1},'fontsize',13)
text(5,-0.3,Dim_names{2,1}{2},'fontsize',13)
text(9,-0.3,Dim_names{2,1}{3},'fontsize',13)
text(13,-0.3,Dim_names{2,1}{4},'fontsize',13)
text(17.2,1,Dim_names{2,1}{1},'fontsize',13,'Rotation',-90)
text(17.2,5,Dim_names{2,1}{2},'fontsize',13,'Rotation',-90)
text(17.2,9,Dim_names{2,1}{3},'fontsize',13,'Rotation',-90)
text(17.2,13,Dim_names{2,1}{4},'fontsize',13,'Rotation',-90)
set(h1,'Units','normalized')




r = 2;
c = 2;
tps = [-0.05 0.05 0.1 0.2 0.4 0.6];

figure('Pos', [222 82 1391 898]);
subplot(4,6,1:18)
plot([]); hold on
dat1 = mean(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,1)),3),1);
dat2 = mean(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,2)),3),1);
SEM1 = nanstd(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,1)),3),0,1)./sqrt(size(RSA_Data_Enc.(ROI{r}).red16_Data,1));
SEM2 = nanstd(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,2)),3),0,1)./sqrt(size(RSA_Data_Enc.(ROI{r}).red16_Data,1));
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat1 fliplr(dat1 + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat1 fliplr(dat1 - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat2 fliplr(dat2 + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat2 fliplr(dat2 - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
h1 = plot(TimeVec_Enc, dat1,'b','linewidth',2);
h2 = plot(TimeVec_Enc, dat2,'r','linewidth',2);
hold off
ylabel(RSA_Data_Enc.meas16{msr}); title(ROI_names{r})
xlim([-0.1 0.8]);lg = legend([h1 h2], {[Cat_names{c},' ',Comp_names{c,cp}{1}],[Cat_names{c},' ',Comp_names{c,cp}{2}]}); legend boxoff; set(lg,'FontSize',14)
box off; %ylim([0.42 0.6])
set(gca,'linewidth',2.5,'FontSize',14)
for sbp = 1:length(tps)
    h1 = subplot(4,6,18+sbp,'Units','pixels','xlim',[0 18],'ylim',[0 18]);
    cur_data = zeros(16);
    cur_data(UTM_idx) = squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,dsearchn(RSA_Data_Enc.TimeVec',tps(sbp)),:),1));
    if(clim_wi)
        colb_lim = [0.45 0.57];
    end
    imagesc(cur_data);
    axis square
    grid on
    set(get(gca,'Yruler'),'Minortick',Dim_names{2,2}(1,:))
    set(get(gca,'Xruler'),'Minortick',Dim_names{2,2}(1,:))
    set(gca,'clim',colb_lim,'xtick',Dim_names{1,2}(1,:),'xticklabels',[],...
            'ytick',Dim_names{1,2}(2,:),'yticklabels',[], 'TickLength',[0 0],'XMinorgrid','on','YMinorgrid','on',...
            'gridcolor','w','gridalpha',.9,'minorgridlinestyle','--','minorgridalpha',.5,'MinorGridColor','w')
    xlabel(sprintf('%2.0f ms',tps(sbp)*1000))
end
%pos = get(h1,'pos');
%set(h1,'pos',pos)
text(3,-1.4,Dim_names{1,1}{1}(1:3),'fontsize',9)
text(11,-1.4,Dim_names{1,1}{2}(1:3),'fontsize',9)
text(18.2,3,Dim_names{1,1}{1}(1:3),'fontsize',9,'Rotation',-90)
text(18.2,11,Dim_names{1,1}{2}(1:3),'fontsize',9,'Rotation',-90)
text(1,-0.3,Dim_names{2,1}{1}(1:3),'fontsize',9)
text(5,-0.3,Dim_names{2,1}{2}(1:3),'fontsize',9)
text(9,-0.3,Dim_names{2,1}{3}(1:3),'fontsize',9)
text(13,-0.3,Dim_names{2,1}{4}(1:3),'fontsize',9)
text(17.2,1,Dim_names{2,1}{1}(1:3),'fontsize',9,'Rotation',-90)
text(17.2,5,Dim_names{2,1}{2}(1:3),'fontsize',9,'Rotation',-90)
text(17.2,9,Dim_names{2,1}{3}(1:3),'fontsize',9,'Rotation',-90)
text(17.2,13,Dim_names{2,1}{4}(1:3),'fontsize',9,'Rotation',-90)
%set(h1,'Units','normalized')

%saveas(gcf,sprintf('Results/Enc_%s_16_%s_LDA_%s.png',ROI{r},Cat{c},datanames{d}))
%close(gcf)



%% Plot Electrode Activations


r = 4;
elec_idx = ROI_idx{r};

Sensor_Data = zeros(size(RSA_Data_Enc.(ROI{r}).red16_SensorData,1), length(ROI_all_idx), size(RSA_Data_Enc.(ROI{r}).red16_SensorData,3),4);
Sensor_Data(:,elec_idx,:,1) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,1),4);
Sensor_Data(:,elec_idx,:,2) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,2),4);
Sensor_Data(:,elec_idx,:,3) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,3),4);
Sensor_Data(:,elec_idx,:,4) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,4),4);

GA_Sensor_Data = [];
GA_Sensor_Data.avg = squeeze(mean(Sensor_Data(:,:,:,2),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data.var = squeeze(std(Sensor_Data(:,:,:,2),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data.dof = size(Sensor_Data,1).*ones(size(Sensor_Data,2),size(Sensor_Data,3)); 
GA_Sensor_Data.time = RSA_Data_Enc.TimeVec; 
GA_Sensor_Data.label = elecs; 
GA_Sensor_Data.dimord = 'chan_time';



cfg = [];
cfg.xlim             = [0.11 0.17]; %[0.24:0.02:0.54]; %[0.06:0.02:0.36];
cfg.zlim             = prctile(GA_Sensor_Data.avg(:),[2.5 97.5]);
%cfg.colorbar         = 'yes';
cfg.layout           = 'biosemi128.lay';
cfg.highlight        = 'on';
cfg.highlightchannel = elec_idx;
cfg.highlightsymbol  = 's';
cfg.highlightcolor   = [0 0 0];
cfg.comment          = 'xlim';
cfg.commentpos       = 'title';


%figure('pos',[21 105 884 843])
figure
ft_topoplotER(cfg, GA_Sensor_Data);


r = 2;
elec_idx = ROI_idx{r};

Sensor_Data = zeros(size(RSA_Data_Enc.(ROI{r}).red16_SensorData,1), length(ROI_all_idx), size(RSA_Data_Enc.(ROI{r}).red16_SensorData,3),4);
Sensor_Data(:,elec_idx,:,1) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,1),4);
Sensor_Data(:,elec_idx,:,2) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,2),4);
Sensor_Data(:,elec_idx,:,3) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,3),4);
Sensor_Data(:,elec_idx,:,4) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,4),4);


GA_Sensor_Data = [];
GA_Sensor_Data.avg = squeeze(mean(Sensor_Data(:,:,:,4),1)); %squeeze(mean(Sensor_Data(:,:,:,4) - Sensor_Data(:,:,:,3),1)); 
GA_Sensor_Data.var = squeeze(std(Sensor_Data(:,:,:,4),0,1)); %squeeze(std(Sensor_Data(:,:,:,4) - Sensor_Data(:,:,:,3),0,1));  
GA_Sensor_Data.dof = size(Sensor_Data,1).*ones(size(Sensor_Data,2),size(Sensor_Data,3)); 
GA_Sensor_Data.time = RSA_Data_Enc.TimeVec; 
GA_Sensor_Data.label = elecs; 
GA_Sensor_Data.dimord = 'chan_time';

cfg.xlim             = [0.75 1.2];
cfg.zlim             = prctile(GA_Sensor_Data.avg(:),[2.5 97.5]);
cfg.highlightchannel = elec_idx;

%figure('pos',[916 105 884 843])
figure
ft_topoplotER(cfg, GA_Sensor_Data);





%% Plot RSA Time series
r = 4;
c = 1;
cp = 1;

timelim = [-0.1  1.2];

%TimeFrameTopo = dsearchn(TimeVec_Enc',linspace(timelim(1),timelim(2),8)')';
%TimeFrameTopo = [dsearchn(TimeVec_Enc',-.1) dsearchn(TimeVec_Enc',0) round(linspace(find(sign_mcc_clust > 0,1,'first'),find(sign_mcc_clust > 0,1,'last'),6))];
%TimeFrameTopo = dsearchn(TimeVec_Enc',[0.05 0.08 0.12 0.13 0.14 0.15 0.28 0.4]')'; % Parietal Perceptual Timepoints
%TimeFrameTopo = dsearchn(TimeVec_Enc',[0.12 0.25 0.28 0.37 0.425 0.46 0.57 0.65]')'; % Parietal Semantic Timepoints
%TimeFrameTopo = dsearchn(TimeVec_Enc',[0.05 0.08 0.12 0.13 0.14 0.15 0.28 0.4]')'; % ALL Perceptual Timepoints
TimeFrameTopo = dsearchn(TimeVec_Enc',[0.12 0.28 0.37 0.46 0.7 0.83 0.91 1.1]')'; % ALL Semantic Timepoints

elec_idx = ROI_idx{r};

Sensor_Data = zeros(size(RSA_Data_Enc.(ROI{r}).red16_SensorData,1), length(ROI_all_idx), size(RSA_Data_Enc.(ROI{r}).red16_SensorData,3),4);
Sensor_Data(:,elec_idx,:,1) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,1),4);
Sensor_Data(:,elec_idx,:,2) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,2),4);
Sensor_Data(:,elec_idx,:,3) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,3),4);
Sensor_Data(:,elec_idx,:,4) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,4),4);

GA_Sensor_Data = [];
GA_Sensor_Data.avg1_1 = squeeze(mean(Sensor_Data(:,:,:,1),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data.var1_1 = squeeze(std(Sensor_Data(:,:,:,1),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data.avg1_2 = squeeze(mean(Sensor_Data(:,:,:,2),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data.var1_2 = squeeze(std(Sensor_Data(:,:,:,2),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1)); 
GA_Sensor_Data.avg2_1 = squeeze(mean(Sensor_Data(:,:,:,3),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data.var2_1 = squeeze(std(Sensor_Data(:,:,:,3),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data.avg2_2 = squeeze(mean(Sensor_Data(:,:,:,4),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data.var2_2 = squeeze(std(Sensor_Data(:,:,:,4),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data.dof = size(Sensor_Data,1).*ones(size(Sensor_Data,2),size(Sensor_Data,3)); 
GA_Sensor_Data.time = RSA_Data_Enc.TimeVec; 
GA_Sensor_Data.label = elecs; 
GA_Sensor_Data.dimord = 'chan_time';

RDM_Data = squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data,1));
colb_lim = prctile(RDM_Data(:),[2.5 97.5]);
Dim_names = {};
[a,b] = unique(RSA_Data_Enc.TrialInfo{1,1}(1:64,6),'stable');
Dim_names{1,1} = RSA_Data_Enc.TrialInfo{1,1}([b;b+64],6);
Dim_names{1,2} = [8.5; 8.5];
[a,b] = unique(RSA_Data_Enc.TrialInfo{1,1}(1:64,8),'stable');
Dim_names{2,1} = RSA_Data_Enc.TrialInfo{1,1}([b;b+64],8);
Dim_names{2,2} = [4.5:4:12.5; 4.5:4:12.5];

Hyp_idx{1}  = RSA_Data_Enc.Cond_idx.Perc_WI_BT_idx;
Hyp_idx{2}    = RSA_Data_Enc.Cond_idx.Sem_WI_BT_idx;
cfg = [];
cfg.nPerms = 1000;
cfg.thresh_pval = 0.05;
cfg.mcc_cluster_pval = 0.05;
cfg.TimeVec = RSA_Data_Enc.TimeVec;
cfg.Hyp_idx = Hyp_idx{c};
cfg.matshuff = false;
cfg.twoside = true;
Results1 = rsa_perm(cfg, RSA_Data_Enc.(ROI{r}).red16_Data);

cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.center = 'no';
cfg.projection  = 'polar';
cfg.overlap     = 'shift';
cfg.skipscale   = 'yes';
cfg.skipcomnt   = 'yes';
layout_elec = ft_prepare_layout(cfg);

figpos    = [330  52  1256  939];
tspos     = [(1/5)*figpos(3) (2/3)*figpos(4) (3/5)*figpos(3) 280];
padd = 5;
topobtpos = [padd   tspos(2)-figpos(3)/8-60       figpos(3)/8-2*padd    figpos(3)/8-2*padd];
topowipos = [padd   tspos(2)-2*(figpos(3)/8+30)   figpos(3)/8-2*padd    figpos(3)/8-2*padd];
rdmpos    = [padd   tspos(2)-3*(figpos(3)/8+40)   figpos(3)/8-2*padd    figpos(3)/8-2*padd];
roipos    = [tspos(1)+tspos(3)+40  tspos(2)+50  figpos(3)-tspos(1)-tspos(3)-80 figpos(3)-tspos(1)-tspos(3)-80];


figure('Pos', figpos);
tsh = axes('units','pixels','Position',tspos);
hold on
dat1 = mean(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,1)),3),1);
dat2 = mean(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,2)),3),1);
SEM1 = nanstd(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,1)),3),0,1)./sqrt(size(RSA_Data_Enc.(ROI{r}).red16_Data,1));
SEM2 = nanstd(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,2)),3),0,1)./sqrt(size(RSA_Data_Enc.(ROI{r}).red16_Data,1));
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat1 fliplr(dat1 + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat1 fliplr(dat1 - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat2 fliplr(dat2 + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat2 fliplr(dat2 - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
h1 = plot(TimeVec_Enc, dat1,'b','linewidth',2);
h2 = plot(TimeVec_Enc, dat2,'r','linewidth',2);
xlim(timelim)
if(Results1.H == 1)
    sign_mcc_clust = abs(Results1.zmapthresh);
    sign_mcc_clust(sign_mcc_clust > 0) = min(get(gca,'ylim'))*1.02;
    plot(TimeVec_Enc,sign_mcc_clust,'ko','MarkerFaceColor','k')
end
hold off
ylabel(RSA_Data_Enc.meas16{msr}); title([ROI{r}])
lg = legend([h1 h2], {[Cat_names{c},' ',Comp_names{c,cp}{1}],[Cat_names{c},' ',Comp_names{c,cp}{2}]}); legend boxoff; set(lg,'FontSize',14) %,'Location','east')
box off; %ylim([0.42 0.6])
set(gca,'linewidth',2.5,'FontSize',14,'xtick',-0.2:0.1:1.5)

xlim_idx = dsearchn(TimeVec_Enc',tsh.XLim')';
pixpertp = (tsh.Position(3)/diff(xlim_idx));


roi_h = axes('Units','Pixel','Position',roipos);
ft_plot_layout(layout_elec,'label','no','box','no','pointsymbol','o','pointcolor','k','pointsize',3)
h = get(gca,'Child');
hold on
plot(h(6).XData(elec_idx), h(6).YData(elec_idx),'s','MarkerSize',8,'MarkerFaceColor','k')
hold off

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
%cfg.commentpos       = 'middlebottom';

for to = 1:length(TimeFrameTopo)
    topobt_h(to) = axes('units','pixels','Position',topobtpos + (to-1)*[topobtpos(3)+2*padd 0 0 0]);
    cfg.parameter = ['avg',num2str(c),'_2'];
    cfg.xlim      = [TimeVec_Enc(TimeFrameTopo(to)) TimeVec_Enc(TimeFrameTopo(to))];
    cfg.zlim      = prctile(GA_Sensor_Data.(['avg',num2str(c),'_',num2str(cp)])(:),[2.5 97.5]);
    ft_topoplotER(cfg, GA_Sensor_Data);
    h = get(gca,'Children');
    h(1).Color = [1 0 0]; h(2).Color = [1 0 0]; h(3).Color = [1 0 0]; h(4).Color = [1 0 0]; h(5).Color = [1 0 0];
    topowi_h(to) = axes('units','pixels','Position',topowipos + (to-1)*[topowipos(3)+2*padd 0 0 0]);
    cfg.parameter = ['avg',num2str(c),'_1'];
    cfg.xlim      = [TimeVec_Enc(TimeFrameTopo(to)) TimeVec_Enc(TimeFrameTopo(to))];
    cfg.zlim      = prctile(GA_Sensor_Data.(['avg',num2str(c),'_',num2str(cp)])(:),[2.5 97.5]);
    ft_topoplotER(cfg, GA_Sensor_Data);
    h = get(gca,'Children');
    h(1).Color = [0 0 1]; h(2).Color = [0 0 1]; h(3).Color = [0 0 1]; h(4).Color = [0 0 1]; h(5).Color = [0 0 1];
    text(0, min(get(gca,'ylim'))*1.15, [num2str(round(TimeVec_Enc(TimeFrameTopo(to))*1000)),' ms'], 'Fontsize', 15, 'HorizontalAlignment','center')
    rdm_h(to) = axes('units','pixels','Position',rdmpos + (to-1)*[rdmpos(3)+2*padd 0 0 0]);
    cur_data = nan(16);
    cur_data(UTM_idx) = RDM_Data(TimeFrameTopo(to),:);
    imh = imagesc(cur_data);
    axis('off', 'image')
    set(gca,'clim',colb_lim,'xtick',Dim_names{1,2}(1,:),'xticklabels',[],...
            'ytick',Dim_names{1,2}(2,:),'yticklabels',[], 'TickLength',[0 0])
    set(imh, 'AlphaData', ~isnan(cur_data))
    hold on
    if(c == 1)
        triplot([1 2 3],[0.55;8.45;8.45],[0.55;8.45;0.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[0.55 0.55 7.9 7.9],'EdgeColor',[0 0 1],'LineWidth',1)
        triplot([1 2 3],[8.55;16.45;16.45],[8.55;16.45;8.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[8.55 8.55 7.9 7.9],'EdgeColor',[0 0 1],'LineWidth',1)
        rectangle('Position',[8.55 0.55 7.9 7.9],'EdgeColor',[1 0 0],'LineWidth',1)
    else
        triplot([1 2 3],[0.55;4.45;4.45],[0.55;4.45;0.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[ 0.55  0.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
        rectangle('Position',[ 4.55  0.55 3.9 3.9],'EdgeColor',[1 0 0],'LineWidth',1)
        rectangle('Position',[ 8.55  0.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
        rectangle('Position',[12.55  0.55 3.9 3.9],'EdgeColor',[1 0 0],'LineWidth',1)
        triplot([1 2 3],[4.55;8.45;8.45],[4.55;8.45;4.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[ 4.55  4.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
        rectangle('Position',[ 8.55  4.55 3.9 3.9],'EdgeColor',[1 0 0],'LineWidth',1)
        rectangle('Position',[12.55  4.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
        triplot([1 2 3],[8.55;12.45;12.45],[8.55;12.45;8.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[ 8.55  8.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
        rectangle('Position',[12.55  8.55 3.9 3.9],'EdgeColor',[1 0 0],'LineWidth',1)
        triplot([1 2 3],[12.55;16.45;16.45],[12.55;16.45;12.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[12.55 12.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
    end
    hold off
    if(to == 1)
        if(c == 1)
            text(4.5,-0.6,Dim_names{1,1}{1},'fontsize',10, 'HorizontalAlignment','center')
            text(12.5,-0.6,Dim_names{1,1}{2},'fontsize',10, 'HorizontalAlignment','center')
        else
            text(2.5,-0.6,Dim_names{2,1}{1}(1:2),'fontsize',10, 'HorizontalAlignment','center')
            text(6.5,-0.6,Dim_names{2,1}{2}(1:2),'fontsize',10, 'HorizontalAlignment','center')
            text(10.5,-0.6,Dim_names{2,1}{3}(1:2),'fontsize',10, 'HorizontalAlignment','center')
            text(14.5,-0.6,Dim_names{2,1}{4}(1:2),'fontsize',10, 'HorizontalAlignment','center')
        end
    end
    
    annotation('line',(([TimeFrameTopo(to) TimeFrameTopo(to)]-xlim_idx(1))*pixpertp + tsh.Position(1))./figpos(3),[tsh.Position(2)-20 tsh.Position(2)]./figpos(4))
    annotation('line',(topobt_h(to).Position(1)+[topobt_h(to).Position(3)  topobt_h(to).Position(3)]/2)./figpos(3), [topobt_h(to).Position(2)+topobt_h(to).Position(4) topobt_h(to).Position(2)+topobt_h(to).Position(4)+20]./figpos(4))
    annotation('line', [topobt_h(to).Position(1)+topobt_h(to).Position(3)/2  (TimeFrameTopo(to)-xlim_idx(1))*pixpertp + tsh.Position(1)]./figpos(3), [topobt_h(to).Position(2)+topobt_h(to).Position(4)+20 tsh.Position(2)-20]./figpos(4))
    
end



%% Plot for SFN-Poster

r = 4;
c = 2;
cp = 1;

timelim = [-0.1  1.2];

%TimeFrameTopo = dsearchn(TimeVec_Enc',linspace(timelim(1),timelim(2),8)')';
%TimeFrameTopo = [dsearchn(TimeVec_Enc',-.1) dsearchn(TimeVec_Enc',0) round(linspace(find(sign_mcc_clust > 0,1,'first'),find(sign_mcc_clust > 0,1,'last'),6))];
%TimeFrameTopo = dsearchn(TimeVec_Enc',[0.05 0.08 0.12 0.13 0.14 0.15 0.28 0.4]')'; % Parietal Perceptual Timepoints
%TimeFrameTopo = dsearchn(TimeVec_Enc',[0.12 0.25 0.28 0.37 0.425 0.46 0.57 0.65]')'; % Parietal Semantic Timepoints
%TimeFrameTopo = dsearchn(TimeVec_Enc',[0.05 0.08 0.12 0.13 0.14 0.15 0.28 0.4]')'; % ALL Perceptual Timepoints
TimeFrameTopo = dsearchn(TimeVec_Enc',[0.12 0.28 0.37 0.46 0.7 0.83 0.91 1.1]')'; % ALL Semantic Timepoints

elec_idx = ROI_idx{r};

Sensor_Data = zeros(size(RSA_Data_Enc.(ROI{r}).red16_SensorData,1), length(ROI_all_idx), size(RSA_Data_Enc.(ROI{r}).red16_SensorData,3),4);
Sensor_Data(:,elec_idx,:,1) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,1),4);
Sensor_Data(:,elec_idx,:,2) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,2),4);
Sensor_Data(:,elec_idx,:,3) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,3),4);
Sensor_Data(:,elec_idx,:,4) = mean(RSA_Data_Enc.(ROI{r}).red16_SensorData(:,:,:,4),4);

GA_Sensor_Data = [];
GA_Sensor_Data.avg1_1 = squeeze(mean(Sensor_Data(:,:,:,1),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data.var1_1 = squeeze(std(Sensor_Data(:,:,:,1),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data.avg1_2 = squeeze(mean(Sensor_Data(:,:,:,2),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data.var1_2 = squeeze(std(Sensor_Data(:,:,:,2),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1)); 
GA_Sensor_Data.avg2_1 = squeeze(mean(Sensor_Data(:,:,:,3),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data.var2_1 = squeeze(std(Sensor_Data(:,:,:,3),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data.avg2_2 = squeeze(mean(Sensor_Data(:,:,:,4),1)); %squeeze(mean(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),1)); 
GA_Sensor_Data.var2_2 = squeeze(std(Sensor_Data(:,:,:,4),0,1)); %squeeze(std(Sensor_Data(:,:,:,2) - Sensor_Data(:,:,:,1),0,1));  
GA_Sensor_Data.dof = size(Sensor_Data,1).*ones(size(Sensor_Data,2),size(Sensor_Data,3)); 
GA_Sensor_Data.time = RSA_Data_Enc.TimeVec; 
GA_Sensor_Data.label = elecs; 
GA_Sensor_Data.dimord = 'chan_time';

RDM_Data = squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data,1));
colb_lim = prctile(RDM_Data(:),[2.5 97.5]);
Dim_names = {};
[a,b] = unique(RSA_Data_Enc.TrialInfo{1,1}(1:64,6),'stable');
Dim_names{1,1} = RSA_Data_Enc.TrialInfo{1,1}([b;b+64],6);
Dim_names{1,2} = [8.5; 8.5];
[a,b] = unique(RSA_Data_Enc.TrialInfo{1,1}(1:64,8),'stable');
Dim_names{2,1} = RSA_Data_Enc.TrialInfo{1,1}([b;b+64],8);
Dim_names{2,2} = [4.5:4:12.5; 4.5:4:12.5];

Hyp_idx{1}  = RSA_Data_Enc.Cond_idx.Perc_WI_BT_idx;
Hyp_idx{2}    = RSA_Data_Enc.Cond_idx.Sem_WI_BT_idx;
cfg = [];
cfg.nPerms = 1000;
cfg.thresh_pval = 0.05;
cfg.mcc_cluster_pval = 0.05;
cfg.TimeVec = RSA_Data_Enc.TimeVec;
cfg.Hyp_idx = Hyp_idx{c};
cfg.matshuff = false;
cfg.twoside = true;
Results1 = rsa_perm(cfg, RSA_Data_Enc.(ROI{r}).red16_Data);

figure('Pos', [300 200 1000 500]);
tsh = axes('units','pixels');
hold on
dat1 = mean(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,1)),3),1);
dat2 = mean(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,2)),3),1);
SEM1 = nanstd(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,1)),3),0,1)./sqrt(size(RSA_Data_Enc.(ROI{r}).red16_Data,1));
SEM2 = nanstd(mean(RSA_Data_Enc.(ROI{r}).red16_Data(:,:,RSA_Data_Enc.Cond_idx.([Cat{c},'_',Comp{cp},'_idx'])(:,2)),3),0,1)./sqrt(size(RSA_Data_Enc.(ROI{r}).red16_Data,1));
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat1 fliplr(dat1 + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat1 fliplr(dat1 - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat2 fliplr(dat2 + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec_Enc fliplr(TimeVec_Enc)],[dat2 fliplr(dat2 - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
h1 = plot(TimeVec_Enc, dat1,'b','linewidth',2);
h2 = plot(TimeVec_Enc, dat2,'r','linewidth',2);
xlim(timelim)
if(Results1.H == 1)
    sign_mcc_clust = abs(Results1.zmapthresh);
    sign_mcc_clust(sign_mcc_clust > 0) = min(get(gca,'ylim'))*1.02;
    plot(TimeVec_Enc,sign_mcc_clust,'ko','MarkerFaceColor','k')
end
hold off
ylabel('LDA Accuracy');
lg = legend([h1 h2], {[Cat_names{c},' ',Comp_names{c,cp}{1}],[Cat_names{c},' ',Comp_names{c,cp}{2}]}); legend boxoff; set(lg,'FontSize',14) %,'Location','east')
box off; %ylim([0.42 0.6])
set(gca,'linewidth',2.5,'FontSize',14,'xtick',-0.2:0.1:1.5)


GA_Sensor_Data_avgbtwi = [];
GA_Sensor_Data_avgbtwi.dof = GA_Sensor_Data.dof;
GA_Sensor_Data_avgbtwi.time = GA_Sensor_Data.time;
GA_Sensor_Data_avgbtwi.label = GA_Sensor_Data.label;
GA_Sensor_Data_avgbtwi.dimord = GA_Sensor_Data.dimord;
GA_Sensor_Data_avgbtwi.avg = (GA_Sensor_Data.avg2_2);
GA_Sensor_Data_avgbtwi.var = (GA_Sensor_Data.var2_2);

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
%cfg.commentpos       = 'middlebottom';

cfg.parameter = 'avg';
cfg.xlim      = [0.42 0.48]; %[0.1 0.17]; %[0.18 0.24]; [0.25 0.31];
%cfg.xlim      = [0.460 0.460]; %[0.370 0.370]; %[0.420 0.420]; [0.460 0.460];
cfg.zlim      = prctile(GA_Sensor_Data_avgbtwi.avg(:),[2.5 97.5]);
ft_topoplotER(cfg, GA_Sensor_Data_avgbtwi);


figure
rdm_h(1) = axes('units','pixels');
cur_data = nan(16);
cur_data(UTM_idx) = RDM_Data(250,:);
imh = imagesc(cur_data);
axis('off', 'image')
set(gca,'clim',colb_lim,'xtick',Dim_names{1,2}(1,:),'xticklabels',[],...
        'ytick',Dim_names{1,2}(2,:),'yticklabels',[], 'TickLength',[0 0])
set(imh, 'AlphaData', ~isnan(cur_data))
hold on
if(c == 1)
    triplot([1 2 3],[0.55;8.45;8.45],[0.55;8.45;0.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[0.55 0.55 7.9 7.9],'EdgeColor',[0 0 1],'LineWidth',1)
    triplot([1 2 3],[8.55;16.45;16.45],[8.55;16.45;8.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[8.55 8.55 7.9 7.9],'EdgeColor',[0 0 1],'LineWidth',1)
    rectangle('Position',[8.55 0.55 7.9 7.9],'EdgeColor',[1 0 0],'LineWidth',1)
else
    triplot([1 2 3],[0.55;4.45;4.45],[0.55;4.45;0.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[ 0.55  0.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
    rectangle('Position',[ 4.55  0.55 3.9 3.9],'EdgeColor',[1 0 0],'LineWidth',1)
    rectangle('Position',[ 8.55  0.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
    rectangle('Position',[12.55  0.55 3.9 3.9],'EdgeColor',[1 0 0],'LineWidth',1)
    triplot([1 2 3],[4.55;8.45;8.45],[4.55;8.45;4.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[ 4.55  4.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
    rectangle('Position',[ 8.55  4.55 3.9 3.9],'EdgeColor',[1 0 0],'LineWidth',1)
    rectangle('Position',[12.55  4.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
    triplot([1 2 3],[8.55;12.45;12.45],[8.55;12.45;8.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[ 8.55  8.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
    rectangle('Position',[12.55  8.55 3.9 3.9],'EdgeColor',[1 0 0],'LineWidth',1)
    triplot([1 2 3],[12.55;16.45;16.45],[12.55;16.45;12.55],'Color',[0 0 1],'LineWidth',1) %rectangle('Position',[12.55 12.55 3.9 3.9],'EdgeColor',[0 0 1],'LineWidth',1)
end
hold off
if(c == 1)
    text(4.5,-0.6,Dim_names{1,1}{1},'fontsize',10, 'HorizontalAlignment','center')
    text(12.5,-0.6,Dim_names{1,1}{2},'fontsize',10, 'HorizontalAlignment','center')
else
    text(2.5,-0.6,Dim_names{2,1}{1}(1:2),'fontsize',10, 'HorizontalAlignment','center')
    text(6.5,-0.6,Dim_names{2,1}{2}(1:2),'fontsize',10, 'HorizontalAlignment','center')
    text(10.5,-0.6,Dim_names{2,1}{3}(1:2),'fontsize',10, 'HorizontalAlignment','center')
    text(14.5,-0.6,Dim_names{2,1}{4}(1:2),'fontsize',10, 'HorizontalAlignment','center')
end



figure
rdm_h(1) = axes('units','pixels');
cur_data = nan(16);
cur_data(UTM_idx) = RDM_Data(52,:);
imh = imagesc(cur_data(end-16:end-1));
set(gca,'xtick','','ytick','')




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