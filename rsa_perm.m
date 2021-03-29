%% RSA Permutation Test

function Results = rsa_perm(cfg, Data)

nPerms = cfg.nPerms;
thresh_pval = cfg.thresh_pval;
mcc_cluster_pval = cfg.mcc_cluster_pval;
TimeVec = cfg.TimeVec;
Hyp_Mat = cfg.Hyp_Mat;
matshuff = cfg.matshuff;
twoside = cfg.twoside;

if(twoside)
    ts_os_fac = 0.5;
else
    ts_os_fac = 1;
end

%% Real Data Clusters

grp1 = nanmean(Data(:,:,Hyp_Mat(:) == -1),3);
grp2 = nanmean(Data(:,:,Hyp_Mat(:) == 1),3);
% compute actual t-test of difference (using unequal N and std)
real_diff   = nanmean(grp1 - grp2,1);
tdenom = nanstd(grp1 - grp2,0,1)./sqrt(size(grp1,1));
Real_T = real_diff./tdenom;


%% Surrogate Data Cluster

Surr_T = zeros(nPerms,length(TimeVec));
max_pixel_pvals = zeros(nPerms,2);
max_clust_info_pos = zeros(nPerms,2); max_clust_info_neg = zeros(nPerms,2);
h = waitbar(0,'');
for p = 1:nPerms
    
    if(matshuff)
        % Create a shuffled version of the Hypothesis Matrix
        Shuff_Mat = zeros(size(Hyp_Mat));
        Ind_Mat = triu(reshape(1:size(Hyp_Mat,1)^2,size(Hyp_Mat,1),size(Hyp_Mat,2)));
        rand_idx = randperm(size(Hyp_Mat,1));
        for row = 1:size(Hyp_Mat,1)-1
            for col = (row+1):size(Hyp_Mat,1)
                if(Ind_Mat(rand_idx(row),rand_idx(col)) ~= 0)
                    Shuff_Mat(row,col) = Ind_Mat(rand_idx(row),rand_idx(col));
                else
                    Shuff_Mat(row,col) = Ind_Mat(rand_idx(col),rand_idx(row));
                end
            end
        end

        %curData = permute(Data,[4 3 1 2]);
        grp1 = nanmean(Data(:,:,Shuff_Mat(Hyp_Mat(:) == -1)),3);
        grp2 = nanmean(Data(:,:,Shuff_Mat(Hyp_Mat(:) == 1)),3);
    else
        grp1 = zeros(size(Data,1),size(Data,2)); 
        grp2 = zeros(size(Data,1),size(Data,2));
        conditions = [ones(floor(0.5*size(Data,1)),1); zeros(ceil(0.5*size(Data,1)),1)];
        rand_idx = logical(conditions(randperm(length(conditions))));
        %rand_idx = rand(size(Data,1),1) >= .5;
        grp1( rand_idx,:) = nanmean(Data( rand_idx,:,Hyp_Mat(:) == -1),3);
        grp1(~rand_idx,:) = nanmean(Data(~rand_idx,:,Hyp_Mat(:) ==  1),3);
        grp2( rand_idx,:) = nanmean(Data( rand_idx,:,Hyp_Mat(:) ==  1),3);
        grp2(~rand_idx,:) = nanmean(Data(~rand_idx,:,Hyp_Mat(:) == -1),3);
    end
    % compute actual t-test of difference (using unequal N and std)
    surr_diff   = squeeze(nanmean(grp1 - grp2,1));
    tdenom = nanstd(grp1 - grp2,0,1)./sqrt(size(grp1,1));
    Surr_T(p,:) = surr_diff./tdenom;
    
    
    % save maximum pixel values
    max_pixel_pvals(p,:) = [min(Surr_T(p,:)) max(Surr_T(p,:))];
    
    % Get positive clusters
    pos_clustmap = Surr_T(p,:);
    pos_clustmap(pos_clustmap < tinv(1-thresh_pval*ts_os_fac, size(Data,1)-1)) = 0;
    neg_clustmap = Surr_T(p,:);
    neg_clustmap(neg_clustmap > tinv(thresh_pval*ts_os_fac, size(Data,1)-1)) = 0;

    % get number of elements in largest supra-threshold cluster
    clust_struct_pos = bwconncomp(pos_clustmap);
    if(~isempty(clust_struct_pos.PixelIdxList))
        max_clust_info_pos(p,:) = [max(cellfun(@length, clust_struct_pos.PixelIdxList))   sum(pos_clustmap(clust_struct_pos.PixelIdxList{find(cellfun(@length, clust_struct_pos.PixelIdxList) == max(cellfun(@length, clust_struct_pos.PixelIdxList)),1,'first')}))];
    else
        max_clust_info_pos(p,:) = [0 0];
    end

    clust_struct_neg = bwconncomp(neg_clustmap);
    if(~isempty(clust_struct_neg.PixelIdxList))
        max_clust_info_neg(p,:) = [max(cellfun(@length, clust_struct_neg.PixelIdxList))   sum(neg_clustmap(clust_struct_neg.PixelIdxList{find(cellfun(@length, clust_struct_neg.PixelIdxList) == max(cellfun(@length, clust_struct_neg.PixelIdxList)),1,'first')}))];
    else
        max_clust_info_neg(p,:) = [0 0];
    end
    clear clust_struct_pos clust_struct_neg pos_clustmap neg_clustmap
    
    waitbar(p/nPerms,h,sprintf('%d of %d Permutation finished!',p,nPerms))
end
close(h)


%% Z Map thresholded

% Plot significant uncorrected areas
zmapthresh_pos = Real_T;
zmapthresh_pos(zmapthresh_pos < tinv(1-thresh_pval*ts_os_fac, size(Data,1)-1)) = 0;
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

zmapthresh_neg = Real_T;
zmapthresh_neg(zmapthresh_neg > tinv(thresh_pval*ts_os_fac, size(Data,1)-1)) = 0;
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

zmapthresh = zmapthresh_pos + zmapthresh_neg;
zmapthresh(zmapthresh == 0) = NaN;

%% Save into Results

Results.nPerms  = nPerms;
Results.TimeVec = TimeVec;
Results.Hyp_Mat = Hyp_Mat;
Results.Real_T  = Real_T;
Results.Surr_T  = Surr_T;
Results.max_pixel_pvals  = max_pixel_pvals;
Results.max_clust_info_pos  = max_clust_info_pos;
Results.clustinfo_pos = clustinfo_pos;
Results.max_clust_info_neg  = max_clust_info_neg;
Results.clustinfo_neg = clustinfo_neg;
Results.zmapthresh = zmapthresh;

end