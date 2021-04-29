%% Compute Crosscomparison of RSA timecourse

function CrossComp_RSA = create_xcomp_rsa(cfg, Data1, Data2)

% Parameter Settings
subs = cfg.subs;
slide_window_s = cfg.slide_window; %s
slide_step_s = cfg.slide_step; %s
if(strcmp(cfg.window_average,'gaussian'))
    average_kern = @(x, dim, gsize) sum(bsxfun(@times, x, (gausswin(gsize)./sum(gausswin(gsize)))'),dim);
    cfg.av = 'gaussian'; % Gaussian Average
else
    average_kern = @(x, dim, gsize) mean(x,dim);
    cfg.av = 'uniform';  % Uniform Average
end
timewind1 = cfg.timewind1;
timewind2 = cfg.timewind2;
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
timewind1_idx = dsearchn(TimeVec_Data1',timewind1');
samplingrate1 = unique(1./diff(TimeVec_Data1)); %Hz
slide_window_idx1 = round(samplingrate1*(slide_window_s/2));
slide_step_idx1 = round(samplingrate1*slide_step_s);
TimeVec1_idx = timewind1_idx(1):slide_step_idx1:timewind1_idx(2);
TimeVec1 = TimeVec_Data1(TimeVec1_idx);
cfg.TimeVec1 = TimeVec1;

timewind2_idx = dsearchn(TimeVec_Data2',timewind2');
samplingrate2 = unique(1./diff(TimeVec_Data2)); %Hz
slide_window_idx2 = round(samplingrate2*(slide_window_s/2));
slide_step_idx2 = round(samplingrate2*slide_step_s);
TimeVec2_idx = timewind2_idx(1):slide_step_idx2:timewind2_idx(2);
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


per_ind = Hyp_perceptual{2}(Hyp_perceptual{2}(:) ~= 0);
sem_ind = Hyp_semantic{2}(Hyp_semantic{2}(:) ~= 0);
X1 = kron(eye(length(subs)), [double(per_ind > 0), double(per_ind < 0)]);
X2 = kron(eye(length(subs)), [double(sem_ind > 0), double(sem_ind < 0)]);
Xy = kron(eye(length(subs)), ones(length(per_ind),1));

CrossComp_RSA.RSA_red16   = [];
for r = 1:length(curROI)
    
    Cur_Data1 = Data1.(curROI{r}).red16_Data;
    Cur_Data2 = Data2.(curROI{r}).red16_Data;
    
    if(permtest)
        SurCorr = zeros(n_perms, length(TimeVec2), length(TimeVec1));
        SurMeth1_per = zeros(n_perms, length(TimeVec2), length(TimeVec1));
        SurMeth1_sem = zeros(n_perms, length(TimeVec2), length(TimeVec1));
        
        if(matshuffle)
            %Ind_Mat = triu(reshape(1:size(Hyp_perceptual{2},1)^2,size(Hyp_perceptual{2},1),size(Hyp_perceptual{2},2))).*~eye(16);
            Ind_Mat = triu(ones(size(Hyp_perceptual{2},1))).*~eye(16).*flip(~eye(16));
            Ind_Mat(Ind_Mat(:) > 0) = 1:sum(Ind_Mat(:) > 0);
            All_Shuff_Mat = zeros(112,n_perms);
            for permi = 1:n_perms
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
                All_Shuff_Mat(:,permi) = Shuff_Mat(Shuff_Mat(:) > 0);
            end
        end
    end
    
    %CrossComp_RSA.RSA_red16.(curROI{r}).Corr = zeros(size(Data1.(curROI{r}).red16_Data,1),length(TimeVec2),length(TimeVec1));
    %CrossComp_RSA.RSA_red16.(curROI{r}).Meth1_per = zeros(size(Data1.(curROI{r}).red16_Data,1),length(TimeVec2),length(TimeVec1));
    %CrossComp_RSA.RSA_red16.(curROI{r}).Meth1_sem = zeros(size(Data1.(curROI{r}).red16_Data,1),length(TimeVec2),length(TimeVec1));
    Corr      = zeros(size(Data1.(curROI{r}).red16_Data,1),length(TimeVec2),length(TimeVec1));
    Meth1_per = zeros(size(Data1.(curROI{r}).red16_Data,1),length(TimeVec2),length(TimeVec1));
    Meth1_sem = zeros(size(Data1.(curROI{r}).red16_Data,1),length(TimeVec2),length(TimeVec1));
    
    fprintf('\n')
    nbytes = fprintf('ROI: %s  --  Progress Cross Correlation:  0.0 %%',curROI{r});
    for tp1 = 1:length(TimeVec1)
        
        time_window1 = max(TimeVec1_idx(tp1) - slide_window_idx1, 1):min(TimeVec1_idx(tp1) + slide_window_idx1, length(TimeVec_Data1));
        
        parfor tp2 = 1:length(TimeVec2)
            
            time_window2 = max(TimeVec2_idx(tp2) - slide_window_idx2, 1):min(TimeVec2_idx(tp2) + slide_window_idx2, length(TimeVec_Data2));
            
            curData1 = squeeze(average_kern(Cur_Data1(:,time_window1,Hyp_perceptual{2}(:) ~= 0),2,length(time_window1)))';
            curData2 = squeeze(average_kern(Cur_Data2(:,time_window2,Hyp_perceptual{2}(:) ~= 0),2,length(time_window2)))';
            
            curData1_tr = tiedrank_(curData1,1);
            curData2_tr = tiedrank_(curData2,1);
            
            Corr(:,tp2,tp1) = atanh(fast_corr(curData1_tr,curData2_tr)');
            
            % Weighted Correlation
            %CrossComp_RSA.RSA_red16.(curROI{r}).Meth1_per(:,tp2,tp1) = CrossComp_RSA.RSA_red16.(curROI{r}).Corr(:,tp2,tp1).*weight_corr(curData1_tr, curData2_tr, per_ind);
            %CrossComp_RSA.RSA_red16.(curROI{r}).Meth1_sem(:,tp2,tp1) = CrossComp_RSA.RSA_red16.(curROI{r}).Corr(:,tp2,tp1).*weight_corr(curData1_tr, curData2_tr, sem_ind);
            [W_per, W_sem] = weight_corr(curData1, curData2, X1, X2, Xy);
            Meth1_per(:,tp2,tp1) = Corr(:,tp2,tp1) .* W_per;
            Meth1_sem(:,tp2,tp1) = Corr(:,tp2,tp1) .* W_sem;
            
            if(permtest)
                for permi = 1:n_perms
                    if(matshuffle)
                        surData = curData2(All_Shuff_Mat(:,permi),:); % This makes it consistent shuffling between time points
                        surData_tr = curData2_tr(All_Shuff_Mat(:,permi),:);
                        %surData = tiedrank_(squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,Shuff_Mat(Hyp_perceptual{2}(:) > 0)),2,length(time_window2)))',1);
                    else
                        surData = curData2_tr(randperm(size(curData2_tr,1)),:);
                        %surData = tiedrank_(squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,hyp_idx(randperm(length(hyp_idx)))),2,length(time_window2)))',1);
                    end
                    if(studentized)
                        mu_22 = nanmean((bsxfun(@minus, curData1_tr, nanmean(curData1_tr,1)).^2).*(bsxfun(@minus, surData, nanmean(surData,1)).^2),1);
                        SurCorr(permi,tp2,tp1) = nanmean(atanh(fast_corr(curData1_tr,surData)'./sqrt(bsxfun(@rdivide, mu_22, (nanvar(curData_tr1,1,1).*nanvar(surData,1,1))))));
                    else
                        SurCorr(permi,tp2,tp1) = mean(atanh(fast_corr(curData1_tr,surData_tr)));
                        [W_per, W_sem] = weight_corr(curData1,surData,X1,X2,Xy);
                        SurMeth1_per(permi,tp2,tp1) = mean(atanh(fast_corr(curData1_tr,surData_tr))' .* W_per);
                        SurMeth1_sem(permi,tp2,tp1) = mean(atanh(fast_corr(curData1_tr,surData_tr))' .* W_sem);
                    end
                end
            end
            %fprintf(repmat('\b',1,nbytes))
            %nbytes = fprintf('ROI: %s  --  Progress Cross Correlation: %3.4f %%',curROI{r},((tp2 + (length(TimeVec2))*(tp1-1)) / (length(TimeVec1)*length(TimeVec2)))*100);
            
        end
        fprintf(repmat('\b',1,nbytes))
        nbytes = fprintf('ROI: %s  --  Progress Cross Correlation: %3.2f %%',curROI{r},(tp1 / length(TimeVec1))*100);
    end
    fprintf('\n')
    
    CrossComp_RSA.RSA_red16.(curROI{r}).Corr = Corr;
    CrossComp_RSA.RSA_red16.(curROI{r}).Meth1_per = Meth1_per;
    CrossComp_RSA.RSA_red16.(curROI{r}).Meth1_sem = Meth1_sem;
    
    if(permtest)

        CrossComp_RSA.CorrPermTest.(curROI{r})       = get_sign_cluster(SurCorr, CrossComp_RSA.RSA_red16.(curROI{r}).Corr, thresh_pval, ts_os_fac, mcc_cluster_pval);
        CrossComp_RSA.Meth1_per_PermTest.(curROI{r}) = get_sign_cluster(SurMeth1_per, CrossComp_RSA.RSA_red16.(curROI{r}).Meth1_per, thresh_pval, ts_os_fac, mcc_cluster_pval);
        CrossComp_RSA.Meth1_sem_PermTest.(curROI{r}) = get_sign_cluster(SurMeth1_sem, CrossComp_RSA.RSA_red16.(curROI{r}).Meth1_sem, thresh_pval, ts_os_fac, mcc_cluster_pval);
        
    end
    
end

%% Save Data
CrossComp_RSA.xcomp_dim        = 'sub_rettime_enctime'; %'sub_hyp_rettime_enctime';

end


function [W1, W2] = weight_corr(D1,D2,X1,X2,Xy)

%X1 = kron(eye(size(D1,2)), [double(ind1 > 0), double(ind1 < 0)]);
%X2 = kron(eye(size(D1,2)), [double(ind2 > 0), double(ind2 < 0)]);

% Method 1
%X = [ones(size(D1,1),1), D1(:,sub)];
%y = D2(:,1);
%SS_res1 = y'*(eye(size(D1,1)) - X*((X'*X)\X'))*y;
%X = [ones(size(D1,1),1), D1(:,sub), double(ind > 0)];
%SS_res2 = y'*(eye(size(D1,1)) - X*((X'*X)\X'))*y;
%F = (SS_res1 - SS_res2)/(SS_res2 / (size(D1,1) - 3));
%W(sub) = abs(f2z_bloc(F,1,size(D1,1) - 3))*atanh(fast_corr(D1(:,sub), D2(:,sub)));

% Method 1
y = tiedrank_(D1 .* D2,1);
Y = bsxfun(@times, Xy, y(:));
SS_tot = sum(bsxfun(@minus, y, nanmean(y)).^2,1);
SS_res = sum((Y - X1*((X1'*X1)\X1'*Y)).^2,1);
W1 = atanh(sqrt(max(1 - SS_res./SS_tot,0)))';
SS_res = sum((Y - X2*((X2'*X2)\X2'*Y)).^2,1);
W2 = atanh(sqrt(max(1 - SS_res./SS_tot,0)))';

end


function Results = get_sign_cluster(SurData, RealData, thresh_pval, ts_os_fac, mcc_cluster_pval)

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

Results.real_zdata          = real_zdata;
Results.max_pixel_pvals     = max_pixel_pvals;
Results.max_clust_info_pos  = max_clust_info_pos;
Results.max_clust_info_neg  = max_clust_info_neg;
Results.zmapthresh = zmapthresh_pos + zmapthresh_neg;
Results.zmapthresh(Results.zmapthresh == 0) = NaN;
        
end
