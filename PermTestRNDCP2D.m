%% Random Cutpoints Permutation Test for Crosscompare Data

function Results = PermTestRNDCP2D(cfg, Data)

nPerms = cfg.nPerms;
thresh_pval = cfg.thresh_pval;
mcc_cluster_pval = cfg.mcc_cluster_pval;
if(cfg.twoside)
    ts_os_fac = 0.5;
else
    ts_os_fac = 1;
end

%% Create Surrogate Data by Random Cutpoints

SurData = zeros(nPerms, size(Data,2), size(Data,3));
h = waitbar(0,'');
for permi = 1:nPerms

    RCX = datasample(1:size(Data,3),size(Data,1));
    RCY = datasample(1:size(Data,2),size(Data,1));
    
    tmp = zeros(size(Data));
    for sub = 1:size(Data,1)
        tmp(sub,:,:) = Data(sub,[RCY(sub):end 1:RCY(sub)-1],[RCX(sub):end 1:RCX(sub)-1]);
    end
    
    SurData(permi,:,:) = nanmean(tmp,1);
    
    waitbar(permi/nPerms,h,sprintf('%d of %d Permutation finished!',permi,nPerms))
end


%% Get Surrogate Clusters
surr_zdata = bsxfun(@rdivide, bsxfun(@minus, SurData, nanmean(SurData,1)), nanstd(SurData,0,1));

max_pixel_pvals     = zeros(size(SurData,1), 2);
max_clust_info_pos  = zeros(size(SurData,1), 2);
max_clust_info_neg  = zeros(size(SurData,1), 2);
for permi = 1:nPerms

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


%% Get Real Clusters
real_zdata = squeeze(nanmean(bsxfun(@rdivide, bsxfun(@minus, Data, nanmean(SurData,1)), nanstd(SurData,0,1)),1));
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
if(isempty(~isnan(Results.zmapthresh)))
    Results.H = 1;
else
    Results.H = 0;
end

end
