%% CrossDecoding Encoding and Retrieval Load Data

permtest = true;
n_perms = 100;
thresh_pval = 0.05;
mcc_cluster_pval = 0.05;
ts_os_fac = 0.5;

%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs 
ROI = {'OCC';'TMP';'FRT';'PRT'};
ROI_idx = {ROI_occ_idx; ROI_temp_idx; ROI_front_idx; ROI_pari_idx};

% Subject Names
load('RSA_Data_Enc','Subj_names')

% Gaussian Average across time window and downsample
window_average  = 'gaussian';
slide_window_s  = 0.080; %s
slide_step_s    = 0.010; %s
fwhm            = 0.040; %s
if(strcmp(window_average,'gaussian'))
    average_kern = @(x, dim, gsize, sd) sum(bsxfun(@times, x, (gausswin(gsize,sd)./sum(gausswin(gsize,sd)))'),dim);
    av = 'gaussian'; % Gaussian Average
else
    average_kern = @(x, dim, gsize, sd) mean(x,dim);
    av = 'uniform';  % Uniform Average
end

samplingrate = 1024; %Hz


CrossDec_Enc_Ret = [];
CrossDec_Enc_Ret.window_average = window_average;
CrossDec_Enc_Ret.slide_window = slide_window_s;
CrossDec_Enc_Ret.slide_step = slide_step_s;
CrossDec_Enc_Ret.fwhm = fwhm;
CrossDec_Enc_Ret.ROI = ROI;
CrossDec_Enc_Ret.permtest = permtest;
CrossDec_Enc_Ret.n_perms = n_perms;
CrossDec_Enc_Ret.thresh_pval = thresh_pval;
CrossDec_Enc_Ret.mcc_cluster_pval = mcc_cluster_pval;
CrossDec_Enc_Ret.ts_os_fac = ts_os_fac;
CrossDec_Enc_Ret.TimeX = [];
CrossDec_Enc_Ret.TimeY = [];

fprintf('\n')
nbytes = fprintf('Subject: %s  --  Progress Cross Decoding:  0.0 %%',Subj_names{1});
for sub = 1:length(Subj_names)

    % Load Encoding Data
    load(['Preproc_EEG_Data/Encoding_object_locked/Enc_EEG_Data_',Subj_names{sub}])
    eeg_cfg = [];
    eeg_cfg.Name = Subj_names{sub};
    eeg_cfg.chan_label = elecs;
    eeg_cfg.chan_idx = ROI_all_idx';
    eeg_cfg.Art_corr = false;
    eeg_cfg.BL_corr = 'minBL';
    eeg_cfg.BL_wind = [-0.26 -0.05];
    [preprocEEG_Enc] = preproc_EEGData(eeg_cfg, Enc_Data_EEG);

    % Load Retrieval Data
    load(['Preproc_EEG_Data/Retrieval_response_locked/Ret_EEG_Data_',Subj_names{sub}])
    eeg_cfg = [];
    eeg_cfg.Name = Subj_names{sub};
    eeg_cfg.chan_label = elecs;
    eeg_cfg.chan_idx = ROI_all_idx';
    eeg_cfg.Art_corr = false;
    eeg_cfg.BL_corr = 'demean';
    eeg_cfg.BL_wind = [-4 0.2];
    [preprocEEG_Ret] = preproc_EEGData(eeg_cfg, Ret_Data_EEG);

    % Downsampling
    slide_window_idx = round(samplingrate*(slide_window_s/2));
    slide_step_idx = round(samplingrate*slide_step_s);

    TimeVec_idx = 1:slide_step_idx:length(preprocEEG_Enc.TimeVec1024);
    TimeVecX = preprocEEG_Enc.TimeVec1024(TimeVec_idx);
    preprocEEG_Enc.Data_dwnsmp = zeros(size(preprocEEG_Enc.Data,1), size(preprocEEG_Enc.Data,2), length(TimeVecX));
    for tp = 1:length(TimeVecX) 
        time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, length(preprocEEG_Enc.TimeVec1024));
        preprocEEG_Enc.Data_dwnsmp(:,:,tp) = squeeze(average_kern(permute(preprocEEG_Enc.Data(:,:,time_window),[1 3 2]), 2, length(time_window), length(time_window)/round(samplingrate*(fwhm/2))*.44));
    end

    TimeVec_idx = 1:slide_step_idx:length(preprocEEG_Ret.TimeVec1024);
    TimeVecY = preprocEEG_Ret.TimeVec1024(TimeVec_idx);
    preprocEEG_Ret.Data_dwnsmp = zeros(size(preprocEEG_Ret.Data,1), size(preprocEEG_Ret.Data,2), length(TimeVecY));
    for tp = 1:length(TimeVecY) 
        time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, length(preprocEEG_Ret.TimeVec1024));
        preprocEEG_Ret.Data_dwnsmp(:,:,tp) = squeeze(average_kern(permute(preprocEEG_Ret.Data(:,:,time_window),[1 3 2]), 2, length(time_window), length(time_window)/round(samplingrate*(fwhm/2))*.44));
    end

    % Save TimeVecs
    CrossDec_Enc_Ret.TimeX = TimeVecX;
    CrossDec_Enc_Ret.TimeY = TimeVecY;
    
    for r = 1:size(ROI,1)
    
        % Initialize Data
        if(sub == 1)
            CrossDec_Enc_Ret.(ROI{r}).Data_Perc = zeros(length(Subj_names),length(TimeVecY),length(TimeVecX));
            CrossDec_Enc_Ret.(ROI{r}).Data_Sem  = zeros(length(Subj_names),length(TimeVecY),length(TimeVecX));
            if(permtest)
                PermData.(ROI{r}).SurData_Perc = zeros(length(Subj_names),n_perms,length(TimeVecY),length(TimeVecX));
                PermData.(ROI{r}).SurData_Sem = zeros(length(Subj_names),n_perms,length(TimeVecY),length(TimeVecX));
            end
        end


        %% Crossdecoding

        % Current Data
        CurData_Enc = preprocEEG_Enc.Data_dwnsmp(~isnan(sum(preprocEEG_Enc.Data_dwnsmp(:,:,1),2)),ROI_idx{r,1},:);
        CurData_Ret = preprocEEG_Ret.Data_dwnsmp(~isnan(sum(preprocEEG_Ret.Data_dwnsmp(:,:,1),2)),ROI_idx{r,end},:);

        % Perceptual Dimension
        Perceptual_Cond_Enc = strcmp(preprocEEG_Enc.TrialInfo(cell2mat(preprocEEG_Enc.TrialInfo(:,1)) > 0,6),'Drawing') + 2.*strcmp(preprocEEG_Enc.TrialInfo(cell2mat(preprocEEG_Enc.TrialInfo(:,1)) > 0,6),'Picture');
        Perceptual_Cond_Ret = strcmp(preprocEEG_Ret.TrialInfo(cell2mat(preprocEEG_Ret.TrialInfo(:,1)) > 0,6),'Drawing') + 2.*strcmp(preprocEEG_Ret.TrialInfo(cell2mat(preprocEEG_Ret.TrialInfo(:,1)) > 0,6),'Picture');

        cfg = [];
        cfg.metric      = 'auc';
        cfg.classifier  = 'lda';
        %cfg.cv          = 'holdout';
        %cfg.p           = 0.3;
        cfg.cv          = 'kfold';
        cfg.k           = 5;
        cfg.repeat      = 5;
        cfg.preprocess  = 'average_samples';
        cfg.feedback    = 0;

        perf_perc = zeros(length(TimeVecY),length(TimeVecX),2);
        [perf_perc(:,:,1)] = mv_classify_timextime(cfg, CurData_Enc, Perceptual_Cond_Enc, CurData_Ret, Perceptual_Cond_Ret)';
        [perf_perc(:,:,2)] = mv_classify_timextime(cfg, CurData_Ret, Perceptual_Cond_Ret, CurData_Enc, Perceptual_Cond_Enc);
        CrossDec_Enc_Ret.(ROI{r}).Data_Perc(sub,:,:) = mean(perf_perc,3);


        % Semantic Dimension
        Semantic_Cond_Enc = strcmp(preprocEEG_Enc.TrialInfo(~isnan(sum(CurData_Enc(:,:,1),2)),8),'Animate') + 2.*strcmp(preprocEEG_Enc.TrialInfo(~isnan(sum(CurData_Enc(:,:,1),2)),8),'Inanimate');
        Semantic_Cond_Ret = strcmp(preprocEEG_Ret.TrialInfo(~isnan(sum(CurData_Ret(:,:,1),2)),8),'Animate') + 2.*strcmp(preprocEEG_Ret.TrialInfo(~isnan(sum(CurData_Ret(:,:,1),2)),8),'Inanimate');

        cfg = [];
        cfg.metric      = 'auc';
        cfg.classifier  = 'lda';
        %cfg.cv          = 'holdout';
        %cfg.p           = 0.3;
        cfg.cv          = 'kfold';
        cfg.k           = 5;
        cfg.repeat      = 5;
        cfg.preprocess  = 'average_samples';
        cfg.feedback    = 0;

        perf_sem = zeros(length(TimeVecY),length(TimeVecX),2);
        [perf_sem(:,:,1)] = mv_classify_timextime(cfg, CurData_Enc, Semantic_Cond_Enc, CurData_Ret, Semantic_Cond_Ret)';
        [perf_sem(:,:,2)] = mv_classify_timextime(cfg, CurData_Ret, Semantic_Cond_Ret, CurData_Enc, Semantic_Cond_Enc);
        CrossDec_Enc_Ret.(ROI{r}).Data_Sem(sub,:,:) = mean(perf_sem,3);
        
        
        if(permtest)
            
            for permi = 1:n_perms
                
                Sur_Perceptual_Cond_Ret = Perceptual_Cond_Ret(randperm(length(Perceptual_Cond_Ret)));
                
                cfg = [];
                cfg.metric      = 'auc';
                cfg.classifier  = 'lda';
                %cfg.cv          = 'holdout';
                %cfg.p           = 0.3;
                cfg.cv          = 'kfold';
                cfg.k           = 5;
                cfg.repeat      = 5;
                cfg.preprocess  = 'average_samples';
                cfg.feedback    = 0;

                sur_perf_perc = zeros(length(TimeVecY),length(TimeVecX),2);
                [sur_perf_perc(:,:,1)] = mv_classify_timextime(cfg, CurData_Enc, Perceptual_Cond_Enc, CurData_Ret, Sur_Perceptual_Cond_Ret)';
                [sur_perf_perc(:,:,2)] = mv_classify_timextime(cfg, CurData_Ret, Sur_Perceptual_Cond_Ret, CurData_Enc, Perceptual_Cond_Enc);
                Sur_mean_perf_perc(permi,:,:) = mean(sur_perf_perc,3);
                
                Sur_Semantic_Cond_Ret = Semantic_Cond_Ret(randperm(length(Semantic_Cond_Ret)));
                
                cfg = [];
                cfg.metric      = 'auc';
                cfg.classifier  = 'lda';
                %cfg.cv          = 'holdout';
                %cfg.p           = 0.3;
                cfg.cv          = 'kfold';
                cfg.k           = 5;
                cfg.repeat      = 5;
                cfg.preprocess  = 'average_samples';
                cfg.feedback    = 0;

                sur_perf_sem = zeros(length(TimeVecY),length(TimeVecX),2);
                [sur_perf_sem(:,:,1)] = mv_classify_timextime(cfg, CurData_Enc, Semantic_Cond_Enc, CurData_Ret, Sur_Semantic_Cond_Ret)';
                [sur_perf_sem(:,:,2)] = mv_classify_timextime(cfg, CurData_Ret, Sur_Semantic_Cond_Ret, CurData_Enc, Semantic_Cond_Enc);
                Sur_mean_perf_sem(permi,:,:) = mean(sur_perf_sem,3);
                
            end
            PermData.(ROI{r}).SurData_Perc(sub,:,:,:) = Sur_mean_perf_perc;
            PermData.(ROI{r}).SurData_Sem(sub,:,:,:) = Sur_mean_perf_sem;
        end
        
        fprintf(repmat('\b',1,nbytes))
        nbytes = fprintf('Subject: %s  --  Progress Cross Decoding:  %3.1f %%',Subj_names{sub},(r / size(ROI,1))*100);
    end
    
end
fprintf('\n')

if(permtest)
    Cat = {'Perc','Sem'};
    for r = 1:size(ROI,1)
        
        for c = 1:2
            
            RealData = CrossDec_Enc_Ret.(ROI{r}).(['Data_',Cat{c}]);
            SurData = squeeze(mean(PermData.(ROI{r}).(['SurData_',Cat{c}]),1));

            surr_zdata = bsxfun(@rdivide, bsxfun(@minus, SurData, nanmean(SurData,1)), nanstd(SurData,0,1));
            %figure 
            %subplot(2,4,1); a = nanmean(SurData,1); hist(a(:), 35);
            %subplot(2,4,5); b = nanstd(SurData,0,1); hist(b(:),35);
            %subplot(2,4,2); hist(SurData(1,:),35);
            %subplot(2,4,6); hist(SurData(1,:),35);
            %subplot(2,4,3); c = bsxfun(@minus, SurData, a); hist(c(:),35);
            %subplot(2,4,7); d = bsxfun(@rdivide, SurData, b); hist(d(:),35);
            %subplot(2,4,4); e = bsxfun(@rdivide, c, b); hist(e(:),35);

            %surr_zdata = bsxfun(@rdivide, SurCorr, nanstd(SurCorr,0,1));
            max_pixel_pvals     = zeros(size(SurData,1), 2);
            max_clust_info_pos  = zeros(size(SurData,1), 2);
            max_clust_info_neg  = zeros(size(SurData,1), 2);
            for permi = 1:size(SurData,1)

                % save maximum pixel values
                max_pixel_pvals(permi,:) = [ min(surr_zdata(permi,:)) max(surr_zdata(permi,:)) ];

                pos_clustmap = squeeze(surr_zdata(permi,:,:));
                pos_clustmap(pos_clustmap < norminv(1-thresh_pval/2)) = 0;
                neg_clustmap = squeeze(surr_zdata(permi,:,:));
                neg_clustmap(neg_clustmap > norminv(thresh_pval/2)) = 0;

                % get number of elements in largest supra-threshold cluster
                clust_struct_pos = bwconncomp(pos_clustmap);
                if(~isempty(clust_struct_pos.PixelIdxList))
                    max_clust_info_pos(permi,:) = [max(cellfun(@length, clust_struct_pos.PixelIdxList))   sum(pos_clustmap(clust_struct_pos.PixelIdxList{find(cellfun(@length, clust_struct_pos.PixelIdxList) == max(cellfun(@length, clust_struct_pos.PixelIdxList)),1,'first')}))];
                else
                    max_clust_info_pos(permi,:) = [0 0];
                end

                clust_struct_neg = bwconncomp(neg_clustmap);
                if(~isempty(clust_struct_neg.PixelIdxList))
                    max_clust_info_neg(permi,:) = [max(cellfun(@length, clust_struct_neg.PixelIdxList))   sum(neg_clustmap(clust_struct_neg.PixelIdxList{find(cellfun(@length, clust_struct_neg.PixelIdxList) == max(cellfun(@length, clust_struct_neg.PixelIdxList)),1,'first')}))];
                else
                    max_clust_info_neg(permi,:) = [0 0];
                end
                clear clust_struct_pos clust_struct_neg pos_clustmap neg_clustmap

            end

            real_zdata = squeeze(bsxfun(@rdivide, bsxfun(@minus, nanmean(RealData,1), nanmean(SurData,1)), nanstd(SurData,0,1)));
            zmapthresh_pos = real_zdata;
            zmapthresh_pos(zmapthresh_pos < norminv(1-thresh_pval*ts_os_fac)) = 0;
            [clustmap,Num] = bwlabel(zmapthresh_pos);  clustinfo_pos = zeros(Num,3);
            for cl = 1:Num
                clustinfo_pos(cl,1) = cl;
                clustinfo_pos(cl,2) = sum(clustmap(:) == cl);
                clustinfo_pos(cl,3) = sum(zmapthresh_pos(clustmap(:) == cl));
            end
            clust_threshold = prctile(max_clust_info_pos(:,2),100-(mcc_cluster_pval*ts_os_fac)*100);
            if(Num > 0)
                for i = 1:size(clustinfo_pos)
                    if(clustinfo_pos(i,3) < clust_threshold) 
                        zmapthresh_pos(clustmap == clustinfo_pos(i,1)) = 0;
                    end
                end
            end
            clear cl i Num clustmap

            zmapthresh_neg = real_zdata;
            zmapthresh_neg(zmapthresh_neg > norminv(thresh_pval*ts_os_fac)) = 0;
            [clustmap,Num] = bwlabel(zmapthresh_neg);  clustinfo_neg = zeros(Num,3);
            for cl = 1:Num
                clustinfo_neg(cl,1) = cl;
                clustinfo_neg(cl,2) = sum(clustmap(:) == cl);
                clustinfo_neg(cl,3) = sum(zmapthresh_neg(clustmap(:) == cl));
            end
            clust_threshold = prctile(max_clust_info_neg(:,2),(mcc_cluster_pval*ts_os_fac)*100);
            if(Num > 0)
                for i = 1:size(clustinfo_neg)
                    if(clustinfo_neg(i,3) > clust_threshold) 
                        zmapthresh_neg(clustmap == clustinfo_neg(i,1)) = 0;
                    end
                end
            end
            clear cl i Num clustmap

            CrossDec_Enc_Ret.(ROI{r}).(['Results_',Cat{c}]).real_zdata          = real_zdata;
            CrossDec_Enc_Ret.(ROI{r}).(['Results_',Cat{c}]).max_pixel_pvals     = max_pixel_pvals;
            CrossDec_Enc_Ret.(ROI{r}).(['Results_',Cat{c}]).max_clust_info_pos  = max_clust_info_pos;
            CrossDec_Enc_Ret.(ROI{r}).(['Results_',Cat{c}]).max_clust_info_neg  = max_clust_info_neg;
            CrossDec_Enc_Ret.(ROI{r}).(['Results_',Cat{c}]).zmapthresh = zmapthresh_pos + zmapthresh_neg;
            CrossDec_Enc_Ret.(ROI{r}).(['Results_',Cat{c}]).zmapthresh(CrossDec_Enc_Ret.(ROI{r}).(['Results_',Cat{c}]).zmapthresh == 0) = NaN;
            if(sum(~isnan(CrossDec_Enc_Ret.(ROI{r}).(['Results_',Cat{c}]).zmapthresh(:))) > 0)
                CrossDec_Enc_Ret.(ROI{r}).(['Results_',Cat{c}]).H = 1;
            else
                CrossDec_Enc_Ret.(ROI{r}).(['Results_',Cat{c}]).H = 0;
            end
        
        end
        
    end
    
end


%% Save Data

save('CrossDec_Enc_Ret','CrossDec_Enc_Ret')



%% Plotting Contourf

%% Perceptual Dimension

%enc_lim = [-0.2  1.5];
%ret_lim = [-2.5  0.2];

% Mean Performance of Perceptual
%figure('Pos',[189  137  796  782])
%contourf(TimeVecX, TimeVecY, mean_perf_perc, 40,'linestyle','none'); colorbar
%tw_dat = mean_perf_perc(dsearchn(TimeVecY', ret_lim(1)):dsearchn(TimeVecY', ret_lim(2)), dsearchn(TimeVecX', enc_lim(1)):dsearchn(TimeVecX', enc_lim(2)));
%caxis(prctile(tw_dat(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim); 
%hold on
%plot([0 0],[TimeVecY(1) TimeVecY(end)],'--w','linewidth',2)
%plot([TimeVecX(1) TimeVecX(end)],[0 0],'--w','linewidth',2)
%hold off

% Both Performance Directions of Perceptual
%figure('Pos',[189 137 1531 782])
%subplot(1,2,1)
%contourf(TimeVecX, TimeVecY, perf_perc(:,:,1), 40,'linestyle','none'); colorbar
%tw_dat = perf_perc(dsearchn(TimeVecY', ret_lim(1)):dsearchn(TimeVecY', ret_lim(2)), dsearchn(TimeVecX', enc_lim(1)):dsearchn(TimeVecX', enc_lim(2)),1);
%caxis(prctile(tw_dat(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim); 
%hold on
%plot([0 0],[TimeVecY(1) TimeVecY(end)],'--w','linewidth',2)
%plot([TimeVecX(1) TimeVecX(end)],[0 0],'--w','linewidth',2)
%hold off
%subplot(1,2,2)
%contourf(TimeVecX, TimeVecY, perf_perc(:,:,2), 40,'linestyle','none'); colorbar
%tw_dat = perf_perc(dsearchn(TimeVecY', ret_lim(1)):dsearchn(TimeVecY', ret_lim(2)), dsearchn(TimeVecX', enc_lim(1)):dsearchn(TimeVecX', enc_lim(2)),2);
%caxis(prctile(tw_dat(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim); 
%hold on
%plot([0 0],[TimeVecY(1) TimeVecY(end)],'--w','linewidth',2)
%plot([TimeVecX(1) TimeVecX(end)],[0 0],'--w','linewidth',2)
%hold off



%% Semantic Dimension

%enc_lim = [-0.2  1.5];
%ret_lim = [-2.5  0.2];

% Mean Performance of Semantic
%figure('Pos',[189  137  796  782])
%contourf(TimeVecX, TimeVecY, mean_perf_sem, 40,'linestyle','none'); colorbar
%tw_dat = mean_perf_sem(dsearchn(TimeVecY', ret_lim(1)):dsearchn(TimeVecY', ret_lim(2)), dsearchn(TimeVecX', enc_lim(1)):dsearchn(TimeVecX', enc_lim(2)));
%caxis(prctile(tw_dat(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim); 
%hold on
%plot([0 0],[TimeVecY(1) TimeVecY(end)],'--w','linewidth',2)
%plot([TimeVecX(1) TimeVecX(end)],[0 0],'--w','linewidth',2)
%hold off

% Both Performance Directions of Semantic
%figure('Pos',[189 137 1531 782])
%subplot(1,2,1)
%contourf(TimeVecX, TimeVecY, perf_sem(:,:,1), 40,'linestyle','none'); colorbar
%tw_dat = perf_sem(dsearchn(TimeVecY', ret_lim(1)):dsearchn(TimeVecY', ret_lim(2)), dsearchn(TimeVecX', enc_lim(1)):dsearchn(TimeVecX', enc_lim(2)),1);
%caxis(prctile(tw_dat(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim); 
%hold on
%plot([0 0],[TimeVecY(1) TimeVecY(end)],'--w','linewidth',2)
%plot([TimeVecX(1) TimeVecX(end)],[0 0],'--w','linewidth',2)
%hold off
%subplot(1,2,2)
%contourf(TimeVecX, TimeVecY, perf_sem(:,:,2), 40,'linestyle','none'); colorbar
%tw_dat = perf_sem(dsearchn(TimeVecY', ret_lim(1)):dsearchn(TimeVecY', ret_lim(2)), dsearchn(TimeVecX', enc_lim(1)):dsearchn(TimeVecX', enc_lim(2)),2);
%caxis(prctile(tw_dat(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim); 
%hold on
%plot([0 0],[TimeVecY(1) TimeVecY(end)],'--w','linewidth',2)
%plot([TimeVecX(1) TimeVecX(end)],[0 0],'--w','linewidth',2)
%hold off


%% Compare Perceptual to Semantic Dimension

%enc_lim = [-0.2  1.5];
%ret_lim = [-4    0.2];

%figure('Pos',[189 137 1531 782])
%subplot(1,2,1)
%contourf(TimeVecX, TimeVecY, mean_perf_perc, 40,'linestyle','none'); colorbar
%tw_dat = mean_perf_perc(dsearchn(TimeVecY', ret_lim(1)):dsearchn(TimeVecY', ret_lim(2)), dsearchn(TimeVecX', enc_lim(1)):dsearchn(TimeVecX', enc_lim(2)));
%caxis(prctile(tw_dat(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim); 
%hold on
%plot([0 0],[TimeVecY(1) TimeVecY(end)],'--w','linewidth',2)
%plot([TimeVecX(1) TimeVecX(end)],[0 0],'--w','linewidth',2)
%hold off
%subplot(1,2,2)
%contourf(TimeVecX, TimeVecY, mean_perf_sem, 40,'linestyle','none'); colorbar
%tw_dat = mean_perf_sem(dsearchn(TimeVecY', ret_lim(1)):dsearchn(TimeVecY', ret_lim(2)), dsearchn(TimeVecX', enc_lim(1)):dsearchn(TimeVecX', enc_lim(2)));
%caxis(prctile(tw_dat(:),[1 99])); set(gca,'xlim', enc_lim, 'ylim', ret_lim); 
%hold on
%plot([0 0],[TimeVecY(1) TimeVecY(end)],'--w','linewidth',2)
%plot([TimeVecX(1) TimeVecX(end)],[0 0],'--w','linewidth',2)
%hold off

