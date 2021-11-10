%% Compute Crossdecoding between EEG timecourse

function CrossDecoding = create_xdecode(cfg, Data1, Data2)


% Parameter Settings
subs = cfg.subs;
slide_window_s = cfg.slide_window; %s
slide_step_s = cfg.slide_step; %s
fwhm = cfg.fwhm; %s
if(strcmp(cfg.window_average,'gaussian'))
    average_kern = @(x, dim, gsize, sd) sum(bsxfun(@times, x, (gausswin(gsize,sd)./sum(gausswin(gsize,sd)))'),dim);
    av = 'gaussian'; % Gaussian Average
else
    average_kern = @(x, dim, gsize, sd) mean(x,dim);
    av = 'uniform';  % Uniform Average
end
MNN = cfg.MNN;
timewind1 = cfg.timewind1;
timewind2 = cfg.timewind2;
Hyp_perceptual = cfg.Hyp_perceptual;
Hyp_semantic = cfg.Hyp_semantic;
TimeVec_Data1 = Data1.TimeVec1024;
TimeVec_Data2 = Data2.TimeVec1024;
ROI = cfg.ROI;
ROI_idx = cfg.ROI_idx;
pair = cfg.pair;
msr = cfg.msr;
permtest = cfg.permtest;
n_perms = cfg.n_perms;
thresh_pval = cfg.thresh_pval;
mcc_cluster_pval = cfg.mcc_cluster_pval;
ts_os_fac = cfg.ts_os_fac;
matshuffle = cfg.matshuffle;
samplingrate = 1024; %Hz


CrossDecoding = [];
CrossDecoding.window_average = cfg.window_average;
CrossDecoding.slide_window = cfg.slide_window;
CrossDecoding.slide_step = cfg.slide_step;
CrossDecoding.fwhm = cfg.fwhm;
CrossDecoding.ROI = cfg.ROI;
CrossDecoding.permtest = cfg.permtest;
CrossDecoding.n_perms = cfg.n_perms;
CrossDecoding.thresh_pval = cfg.thresh_pval;
CrossDecoding.mcc_cluster_pval = cfg.mcc_cluster_pval;
CrossDecoding.ts_os_fac = cfg.ts_os_fac;
CrossDecoding.TimeX = [];
CrossDecoding.TimeY = [];


for r = 1:size(ROI)


    if(MNN)
        
        nonnan_ind = find(cell2mat(Data1.TrialInfo(:,1)) ~= 0);
        Sigma = zeros(16,length(ROI_idx{r}),length(ROI_idx{r}));
        for c = 1:16
            tmp = zeros(size(Data1.Data,3),length(ROI_idx{r}),length(ROI_idx{r}));
            for t = 1:size(Data1.Data,3)
                dat = Data1.Data(intersect(nonnan_ind, (c-1)*8+(1:8)),ROI_idx{r},t);
                tmp(t,:,:) = cov1para(dat(~sum(isnan(dat),2),:));
            end
            Sigma(c,:,:) = nanmean(tmp,1);
        end
        tmp = squeeze(mean(Sigma,1));
        sigma_inv = tmp^-0.5;
        Data1.Data_corr = zeros(size(Data1.Data,1), length(ROI_idx{r}), size(Data1.Data,3));
        for t = 1:size(Data1.Data,3)
            Data1.Data_corr(:,:,t) = Data1.Data(:,ROI_idx{r},t)*sigma_inv;
        end
        
        nonnan_ind = find(cell2mat(Data2.TrialInfo(:,1)) ~= 0);
        Sigma = zeros(16,length(ROI_idx{r}),length(ROI_idx{r}));
        for c = 1:16
            tmp = zeros(size(Data2.Data,3),length(ROI_idx{r}),length(ROI_idx{r}));
            for t = 1:size(Data2.Data,3)
                dat = Data2.Data(intersect(nonnan_ind, (c-1)*8+(1:8)),ROI_idx{r},t);
                tmp(t,:,:) = cov1para(dat(~sum(isnan(dat),2),:));
            end
            Sigma(c,:,:) = nanmean(tmp,1);
        end
        tmp = squeeze(mean(Sigma,1));
        sigma_inv = tmp^-0.5;
        Data2.Data_corr = zeros(size(Data2.Data,1), length(ROI_idx{r}), size(Data2.Data,3));
        for t = 1:size(Data2.Data,3)
            Data2.Data_corr(:,:,t) = Data2.Data(:,ROI_idx{r},t)*sigma_inv;
        end
        
    else
        Data1.Data_corr = Data1.Data(:,ROI_idx{r},:);
        Data2.Data_corr = Data2.Data(:,ROI_idx{r},:);
    end

    
    % Downsampling
    slide_window_idx = round(samplingrate*(slide_window_s/2));
    slide_step_idx = round(samplingrate*slide_step_s);

    TimeVec_idx = 1:slide_step_idx:length(Data1.TimeVec1024);
    TimeVecX = Data1.TimeVec1024(TimeVec_idx);
    Data1.Data_dwnsmp = zeros(size(Data1.Data_corr,1), size(Data1.Data_corr,2), length(TimeVecX));
    for tp = 1:length(TimeVecX) 
        time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, length(Data1.TimeVec1024));
        Data1.Data_dwnsmp(:,:,tp) = squeeze(average_kern(permute(Data1.Data_corr(:,:,time_window),[1 3 2]), 2, length(time_window), length(time_window)/round(samplingrate*(fwhm/2))*.44));
    end

    TimeVec_idx = 1:slide_step_idx:length(Data2.TimeVec1024);
    TimeVecY = Data2.TimeVec1024(TimeVec_idx);
    Data2.Data_dwnsmp = zeros(size(Data2.Data_corr,1), size(Data2.Data_corr,2), length(TimeVecY));
    for tp = 1:length(TimeVecY) 
        time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, length(Data2.TimeVec1024));
        Data2.Data_dwnsmp(:,:,tp) = squeeze(average_kern(permute(Data2.Data_corr(:,:,time_window),[1 3 2]), 2, length(time_window), length(time_window)/round(samplingrate*(fwhm/2))*.44));
    end
    

    % Initialize Data
    CrossDecoding.(ROI{r}).Data = zeros(length(TimeVecY),length(TimeVecX),120);
    if(permtest)
        PermData.(ROI{r}).SurData = zeros(n_perms,length(TimeVecY),length(TimeVecX),120);
    end

    %% Crossdecoding
    
    Cat1_idx = (cell2mat(Data1.TrialInfo(cell2mat(Data1.TrialInfo(:,1)) > 0,5))-1)*8 + cell2mat(Data1.TrialInfo(cell2mat(Data1.TrialInfo(:,1)) > 0,9));
    Cat2_idx = (cell2mat(Data2.TrialInfo(cell2mat(Data2.TrialInfo(:,1)) > 0,5))-1)*8 + cell2mat(Data2.TrialInfo(cell2mat(Data2.TrialInfo(:,1)) > 0,9));
    
    % Current Data
    CurData1 = Data1.Data_dwnsmp(cell2mat(Data1.TrialInfo(:,1)) > 0,:,:);
    CurData2 = Data2.Data_dwnsmp(cell2mat(Data2.TrialInfo(:,1)) > 0,:,:);
    
    ct = 1;
    for i = 1:15
        for j = i+1:16
            
            cfg_cd = [];
            cfg_cd.metric      = 'acc';
            cfg_cd.classifier  = 'lda';
            %cfg_cd.cv          = 'holdout';
            %cfg_cd.p           = 0.3;
            cfg_cd.cv          = 'kfold';
            cfg_cd.k           = 4;
            cfg_cd.repeat      = 5;
            %cfg_cd.preprocess  = 'average_samples';
            cfg_cd.feedback    = 0;
            
            perf_perc = zeros(length(TimeVecY),length(TimeVecX),2);
            [perf_perc(:,:,1)] = mv_classify_timextime(cfg_cd, CurData1(Cat1_idx == i | Cat1_idx == j,:,:), [ones(sum(Cat1_idx == i),1); 2*ones(sum(Cat1_idx == j),1)], CurData2(Cat2_idx == i | Cat2_idx == j,:,:), [ones(sum(Cat2_idx == i),1); 2*ones(sum(Cat2_idx == j),1)])';
            [perf_perc(:,:,2)] = mv_classify_timextime(cfg_cd, CurData2(Cat2_idx == i | Cat2_idx == j,:,:), [ones(sum(Cat2_idx == i),1); 2*ones(sum(Cat2_idx == j),1)], CurData1(Cat1_idx == i | Cat1_idx == j,:,:), [ones(sum(Cat1_idx == i),1); 2*ones(sum(Cat1_idx == j),1)]);
            CrossDecoding.(ROI{r}).Data(:,:,ct) = mean(perf_perc,3);
            ct = ct + 1;
            
        end
    end
    
    fprintf('\n')
    fprintf('Create Crossdecoding %s -- Subject %s -- Meas %i  -- ROI %s !', pair, subs, msr,ROI{r})
    fprintf('\n')
    
    if(permtest)

        for permi = 1:n_perms
            
            ct = 1;
            for i = 1:15
                for j = i+1:16
                    
                    Sur_Perceptual_Cond_Ret = [ones(sum(Cat2_idx == i),1); 2*ones(sum(Cat2_idx == j),1)];
                    Sur_Perceptual_Cond_Ret = Sur_Perceptual_Cond_Ret(randperm(length(Sur_Perceptual_Cond_Ret)));
                    
                    cfg = [];
                    cfg.metric      = 'acc';
                    cfg.classifier  = 'lda';
                    %cfg.cv          = 'holdout';
                    %cfg.p           = 0.3;
                    cfg.cv          = 'kfold';
                    cfg.k           = 4;
                    cfg.repeat      = 5;
                    %cfg.preprocess  = 'average_samples';
                    cfg.feedback    = 0;

                    sur_perf_perc = zeros(length(TimeVecY),length(TimeVecX),2);
                    [sur_perf_perc(:,:,1)] = mv_classify_timextime(cfg, CurData1(Cat1_idx == i | Cat1_idx == j,:,:), [ones(sum(Cat1_idx == i),1); 2*ones(sum(Cat1_idx == j),1)], CurData2(Cat2_idx == i | Cat2_idx == j,:,:), Sur_Perceptual_Cond_Ret)';
                    [sur_perf_perc(:,:,2)] = mv_classify_timextime(cfg, CurData2(Cat2_idx == i | Cat2_idx == j,:,:), Sur_Perceptual_Cond_Ret, CurData1(Cat1_idx == i | Cat1_idx == j,:,:), [ones(sum(Cat1_idx == i),1); 2*ones(sum(Cat1_idx == j),1)]);
                    Sur_mean_perf_perc(permi,:,:,ct) = mean(sur_perf_perc,3);
                    ct = ct + 1;
                    
                end
            end

        end
        PermData.(ROI{r}).SurData(:,:,:,:) = Sur_mean_perf_perc;
        
    end

end

% Save TimeVecs
CrossDecoding.TimeX = TimeVecX;
CrossDecoding.TimeY = TimeVecY;


if(permtest)
    Cat = {'Perc','Sem'};
    for r = 1:size(ROI,1)
        
        for c = 1:2
            
            RealData = CrossDecoding.(ROI{r}).(['Data_',Cat{c}]);
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

            CrossDecoding.(ROI{r}).(['Results_',Cat{c}]).real_zdata          = real_zdata;
            CrossDecoding.(ROI{r}).(['Results_',Cat{c}]).max_pixel_pvals     = max_pixel_pvals;
            CrossDecoding.(ROI{r}).(['Results_',Cat{c}]).max_clust_info_pos  = max_clust_info_pos;
            CrossDecoding.(ROI{r}).(['Results_',Cat{c}]).max_clust_info_neg  = max_clust_info_neg;
            CrossDecoding.(ROI{r}).(['Results_',Cat{c}]).zmapthresh = zmapthresh_pos + zmapthresh_neg;
            CrossDecoding.(ROI{r}).(['Results_',Cat{c}]).zmapthresh(CrossDecoding.(ROI{r}).(['Results_',Cat{c}]).zmapthresh == 0) = NaN;
            if(sum(~isnan(CrossDecoding.(ROI{r}).(['Results_',Cat{c}]).zmapthresh(:))) > 0)
                CrossDecoding.(ROI{r}).(['Results_',Cat{c}]).H = 1;
            else
                CrossDecoding.(ROI{r}).(['Results_',Cat{c}]).H = 0;
            end
        
        end
        
    end
    
end


%% Save Data

%save('CrossDec_Enc_Ret','CrossDec_Enc_Ret')


end