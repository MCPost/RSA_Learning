% Parameter Settings
slide_window_s = cfg.slide_window; %s
slide_step_s = cfg.slide_step; %s
if(strcmp(cfg.window_average,'gaussian'))
    average_kern = @(x, dim, gsize) sum(bsxfun(@times, x, (gausswin(gsize)./sum(gausswin(gsize)))'),dim);
    cfg.av = 'gaussian'; % Gaussian Average
else
    average_kern = @(x, dim, gsize) mean(x,dim);
    cfg.av = 'uniform';  % Uniform Average
end
Hyp_perceptual = cfg.Hyp_perceptual;
Hyp_semantic = cfg.Hyp_semantic;
TimeVec_Data1 = Data1.TimeVec;
TimeVec_Data2 = Data2.TimeVec; 
curROI = cfg.ROI;
only16 = cfg.only16;
permtest = cfg.permtest;
n_perms = cfg.n_perms;
thresh_pval = cfg.thresh_pval;
mcc_cluster_pval = cfg.mcc_cluster_pval;
ts_os_fac = cfg.ts_os_fac;
matshuffle = cfg.matshuffle;
studentized = cfg.studentized;

%% Create RSA Matrix
samplingrate1 = unique(1./diff(TimeVec_Data1)); %Hz
slide_window_idx1 = round(samplingrate1*(slide_window_s/2));
slide_step_idx1 = round(samplingrate1*slide_step_s);
TimeVec1_idx = 1:slide_step_idx1:length(TimeVec_Data1);
TimeVec1 = TimeVec_Data1(TimeVec1_idx);
cfg.TimeVec1 = TimeVec1;

samplingrate2 = unique(1./diff(TimeVec_Data2)); %Hz
slide_window_idx2 = round(samplingrate2*(slide_window_s/2));
slide_step_idx2 = round(samplingrate2*slide_step_s);
TimeVec2_idx = 1:slide_step_idx2:length(TimeVec_Data2);
TimeVec2 = TimeVec_Data2(TimeVec2_idx);
cfg.TimeVec2 = TimeVec2;


% Save Data
CrossComp_RSA = cfg;
CrossComp_RSA.RSA_full   = []; 

fprintf('\n')
fprintf('Create Crosscompare RSA Matrix (Slide Window: %1.3f s | Slide Step: %1.3f s)!', slide_window_s, slide_step_s)
fprintf('\n')

if(~only16)
    
%     per_ind_wi = Hyp_perceptual{1}(:) > 0;
%     per_ind_bt = Hyp_perceptual{1}(:) < 0;
%     sem_ind_wi = Hyp_semantic{1}(:) > 0;
%     sem_ind_bt = Hyp_semantic{1}(:) < 0;
    
    for r = 1:length(curROI)
        
        CrossComp_RSA.RSA_full.(curROI{r}) = zeros(size(Data1.(curROI{r}).full_Data,1),length(TimeVec2),length(TimeVec1));
        fprintf('\n')
        nbytes = fprintf('ROI: %s  --  Progress Cross Correlation:  0.0 %%',curROI{r});
        for tp1 = 1:length(TimeVec1)
            
            time_window1 = max(TimeVec1_idx(tp1) - slide_window_idx1, 1):min(TimeVec1_idx(tp1) + slide_window_idx1, length(TimeVec_Data1));
            
            for tp2 = 1:length(TimeVec2)
                
                time_window2 = max(TimeVec2_idx(tp2) - slide_window_idx2, 1):min(TimeVec2_idx(tp2) + slide_window_idx2, length(TimeVec_Data2));
                
                curData1 = tiedrank_(squeeze(average_kern(Data1.(curROI{r}).full_Data(:,time_window1,Hyp_perceptual{1}(:) > 0),2,length(time_window1)))',1);
                curData2 = tiedrank_(squeeze(average_kern(Data2.(curROI{r}).full_Data(:,time_window2,Hyp_perceptual{1}(:) > 0),2,length(time_window1)))',1);
            
                CrossComp_RSA.RSA_full.(curROI{r})(:,tp2,tp1) = fast_corr(curData1,curData2)';
                %CrossComp_RSA.RSA_full.(curROI{r})(:,2,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).full_Data(:,time_window1,Hyp_perceptual{1}(:) < 0),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).full_Data(:,time_window2,Hyp_perceptual{1}(:) < 0),2,length(time_window2)))', 'type', 'spearman'));
                %CrossComp_RSA.RSA_full.(curROI{r})(:,3,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).full_Data(:,time_window1,Hyp_semantic{1}(:) > 0),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).full_Data(:,time_window2,Hyp_semantic{1}(:) > 0),2,length(time_window2)))', 'type', 'spearman'));
                %CrossComp_RSA.RSA_full.(curROI{r})(:,4,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).full_Data(:,time_window1,Hyp_semantic{1}(:) < 0),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).full_Data(:,time_window2,Hyp_semantic{1}(:) < 0),2,length(time_window2)))', 'type', 'spearman'));
            
            end
            fprintf(repmat('\b',1,nbytes))
            nbytes = fprintf('ROI: %s  --  Progress Cross Correlation: %3.2f %%',curROI{r},(tp1 / length(TimeVec1))*100);
        end
        fprintf('\n') 
    end
    
end


% per_ind_wi = Hyp_perceptual{2}(:) > 0;
% per_ind_bt = Hyp_perceptual{2}(:) < 0;
% sem_ind_wi = Hyp_semantic{2}(:) > 0;
% sem_ind_bt = Hyp_semantic{2}(:) < 0;

CrossComp_RSA.RSA_red16   = [];

for r = 1:length(curROI)
    
    if(permtest)
        SurCorr = zeros(n_perms, length(TimeVec2), length(TimeVec1));
    end
    
    CrossComp_RSA.RSA_red16.(curROI{r}) = zeros(size(Data1.(curROI{r}).red16_Data,1),length(TimeVec2),length(TimeVec1));
    fprintf('\n')
    nbytes = fprintf('ROI: %s  --  Progress Cross Correlation:  0.0 %%',curROI{r});
    for tp1 = 1:length(TimeVec1)
        
        time_window1 = max(TimeVec1_idx(tp1) - slide_window_idx1, 1):min(TimeVec1_idx(tp1) + slide_window_idx1, length(TimeVec_Data1));
        
        for tp2 = 1:length(TimeVec2)
            
            time_window2 = max(TimeVec2_idx(tp2) - slide_window_idx2, 1):min(TimeVec2_idx(tp2) + slide_window_idx2, length(TimeVec_Data2));
            
            curData1 = tiedrank_(squeeze(average_kern(Data1.(curROI{r}).red16_Data(:,time_window1,Hyp_perceptual{2}(:) > 0),2,length(time_window1)))',1);
            curData2 = tiedrank_(squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,Hyp_perceptual{2}(:) > 0),2,length(time_window2)))',1);
            
            CrossComp_RSA.RSA_red16.(curROI{r})(:,tp2,tp1) = fast_corr(curData1,curData2)';
            %CrossComp_RSA.RSA_red16.(curROI{r})(:,2,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).red16_Data(:,time_window1,Hyp_perceptual{2}(:) < 0),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,Hyp_perceptual{2}(:) < 0),2,length(time_window2)))', 'type', 'spearman'));
            %CrossComp_RSA.RSA_red16.(curROI{r})(:,3,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).red16_Data(:,time_window1,Hyp_semantic{2}(:) > 0),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,Hyp_semantic{2}(:) > 0),2,length(time_window2)))', 'type', 'spearman'));
            %CrossComp_RSA.RSA_red16.(curROI{r})(:,4,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).red16_Data(:,time_window1,Hyp_semantic{2}(:) < 0),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,Hyp_semantic{2}(:) < 0),2,length(time_window2)))', 'type', 'spearman'));

            %CrossComp_RSA_Mat(tp_ret,tp_enc) = fast_corr(squeeze(Data1.(ROI{1}).red16_Data(sub,tp_enc,mat_ind)), squeeze(Data2.(ROI{1}).red16_Data(sub,tp_ret,mat_ind)));
            %CrossComp_RSA_Mat(1,tp_ret,tp_enc) = corr(squeeze(Data1.(ROI{2}).red16_Data(sub,tp_enc,mat_ind_wi)), squeeze(Data2.(ROI{2}).red16_Data(sub,tp_ret,mat_ind_wi)), 'type', 'spearman');
            %CrossComp_RSA_Mat(2,tp_ret,tp_enc) = corr(squeeze(Data1.(ROI{2}).red16_Data(sub,tp_enc,mat_ind_bt)), squeeze(Data2.(ROI{2}).red16_Data(sub,tp_ret,mat_ind_bt)), 'type', 'spearman');
            
            if(permtest)
%                 hyp_idx = find(Hyp_perceptual{2}(:) > 0);
%                 for permi = 1:n_perms
%                     if(matshuffle)
%                         Shuff_Mat = zeros(size(Hyp_perceptual{2}));
%                         Ind_Mat = triu(reshape(1:size(Hyp_perceptual{2},1)^2,size(Hyp_perceptual{2},1),size(Hyp_perceptual{2},2)));
%                         rand_idx = randperm(size(Hyp_perceptual{2},1));
%                         for row = 1:size(Hyp_perceptual{2},1)-1
%                             for col = (row+1):size(Hyp_perceptual{2},1)
%                                 if(Ind_Mat(rand_idx(row),rand_idx(col)) ~= 0)
%                                     Shuff_Mat(row,col) = Ind_Mat(rand_idx(row),rand_idx(col));
%                                 else
%                                     Shuff_Mat(row,col) = Ind_Mat(rand_idx(col),rand_idx(row));
%                                 end
%                             end
%                         end
%                         surData = tiedrank_(squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,Shuff_Mat(Hyp_perceptual{2}(:) > 0)),2,length(time_window2)))',1);
%                     else
%                         surData = tiedrank_(squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,hyp_idx(randperm(length(hyp_idx)))),2,length(time_window2)))',1);
%                     end
%                     if(studentized)
%                         mu_22 = nanmean((bsxfun(@minus, curData1, nanmean(curData1,1)).^2).*(bsxfun(@minus, surData, nanmean(surData,1)).^2),1);
%                         SurCorr(permi,tp2,tp1) = nanmean(atanh(fast_corr(curData1,surData)'./sqrt(bsxfun(@rdivide, mu_22, (nanvar(curData1,1,1).*nanvar(surData,1,1))))));
%                     else
%                         SurCorr(permi,tp2,tp1) = nanmean(atanh(fast_corr(curData1,surData)));
%                     end
%                 end

                %Ind_Mat = triu(reshape(1:size(Hyp_perceptual{2},1)^2,size(Hyp_perceptual{2},1),size(Hyp_perceptual{2},2))).*~eye(16);
                Ind_Mat = triu(ones(size(Hyp_perceptual{2},1))).*~eye(16).*flip(~eye(16));
                Ind_Mat(Ind_Mat(:) > 0) = 1:sum(Ind_Mat(:) > 0);
                for permi = 1:n_perms
                    if(matshuffle)
                        Shuff_Mat = zeros(size(Hyp_perceptual{2}));
                        rand_idx = randperm(size(Hyp_perceptual{2},1));
                        for row = 1:size(Hyp_perceptual{2},1)-1
                            for col = (row+1):size(Hyp_perceptual{2},1)
                                if(Ind_Mat(rand_idx(row),rand_idx(col)) ~= 0)
                                    Shuff_Mat(row,col) = Ind_Mat(rand_idx(row),rand_idx(col));
                                else
                                    Shuff_Mat(row,col) = Ind_Mat(rand_idx(col),rand_idx(row));
                                end
                            end
                        end
                        surData = curData2(Shuff_Mat(Shuff_Mat(:) > 0),:);
                        %surData = tiedrank_(squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,Shuff_Mat(Hyp_perceptual{2}(:) > 0)),2,length(time_window2)))',1);
                    else
                        surData = curData2(randperm(size(curData2,1)),:);
                        %surData = tiedrank_(squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,hyp_idx(randperm(length(hyp_idx)))),2,length(time_window2)))',1);
                    end
                    if(studentized)
                        mu_22 = nanmean((bsxfun(@minus, curData1, nanmean(curData1,1)).^2).*(bsxfun(@minus, surData, nanmean(surData,1)).^2),1);
                        SurCorr(permi,tp2,tp1) = nanmean(atanh(fast_corr(curData1,surData)'./sqrt(bsxfun(@rdivide, mu_22, (nanvar(curData1,1,1).*nanvar(surData,1,1))))));
                    else
                        SurCorr(permi,tp2,tp1) = nanmean(atanh(fast_corr(curData1,surData)));
                    end
                end
            end
            
        end
        fprintf(repmat('\b',1,nbytes))
        nbytes = fprintf('ROI: %s  --  Progress Cross Correlation: %3.2f %%',curROI{r},(tp1 / length(TimeVec1))*100);
    end
    fprintf('\n')
    
    if(permtest)
        
        surr_zdata = bsxfun(@rdivide, bsxfun(@minus, SurCorr, nanmean(SurCorr,1)), nanstd(SurCorr,0,1));
        %figure 
        %subplot(2,4,1); a = nanmean(SurCorr,1); hist(a(:), 35);
        %subplot(2,4,5); b = nanstd(SurCorr,0,1); hist(b(:),35);
        %subplot(2,4,2); hist(SurCorr(1,:),35);
        %subplot(2,4,6); hist(SurCorr(1,:),35);
        %subplot(2,4,3); c = bsxfun(@minus, SurCorr, a); hist(c(:),35);
        %subplot(2,4,7); d = bsxfun(@rdivide, SurCorr, b); hist(d(:),35);
        %subplot(2,4,4); e = bsxfun(@rdivide, c, b); hist(e(:),35);
        
        %surr_zdata = bsxfun(@rdivide, SurCorr, nanstd(SurCorr,0,1));
        max_pixel_pvals     = zeros(n_perms, 2);
        max_clust_info_pos  = zeros(n_perms, 2);
        max_clust_info_neg  = zeros(n_perms, 2);
        for permi = 1:n_perms
            
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
        
        real_zdata = squeeze(bsxfun(@rdivide, bsxfun(@minus, nanmean(atanh(CrossComp_RSA.RSA_red16.(curROI{r})),1), nanmean(SurCorr,1)), nanstd(SurCorr,0,1)));
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
        
        CrossComp_RSA.real_zdata{r}          = real_zdata;
        CrossComp_RSA.max_pixel_pvals{r}     = max_pixel_pvals;
        CrossComp_RSA.max_clust_info_pos{r}  = max_clust_info_pos;
        CrossComp_RSA.max_clust_info_neg{r}  = max_clust_info_neg;
        CrossComp_RSA.zmapthresh{r} = zmapthresh_pos + zmapthresh_neg;
        CrossComp_RSA.zmapthresh{r}(CrossComp_RSA.zmapthresh{r} == 0) = NaN;
        
    end
    
end

%% Save Data
CrossComp_RSA.xcomp_dim        = 'sub_rettime_enctime'; %'sub_hyp_rettime_enctime';