%% RSA Permutation Test

function Results = rsa_perm(cfg, Data)

nPerms = cfg.nPerms;
thresh_pval = cfg.thresh_pval;
mcc_cluster_pval = cfg.mcc_cluster_pval;
TimeVec = cfg.TimeVec;
Hyp_Mat = cfg.Hyp_Mat;
matshuff = cfg.matshuff;

%% Real Data Clusters

grp1 = nanmean(Data(:,:,Hyp_Mat(:) == -1),3);
grp2 = nanmean(Data(:,:,Hyp_Mat(:) == 1),3);
% compute actual t-test of difference (using unequal N and std)
real_diff   = nanmean(grp1 - grp2,1);
tdenom = nanstd(grp1 - grp2,0,1)./sqrt(size(grp1,1));
Real_T = real_diff./tdenom;


%% Surrogate Data Cluster

Surr_T = zeros(nPerms,length(TimeVec));
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
        rand_idx = rand(5,1) >= .5;
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
    max_pixel_pvals(p,:) = max(Surr_T(p,:));
    
    % Get positive clusters
    pos_clustmap = Surr_T(p,:);
    pos_clustmap(pos_clustmap < tinv(1-thresh_pval, size(Data,1)-1)) = 0;

    % get number of elements in largest supra-threshold cluster
    [clustmap,Num] = bwlabel(pos_clustmap); clustinfo = zeros(Num,3);
    if(Num > 0)
        for cl = 1:Num
            clustinfo(cl,1) = cl;
            clustinfo(cl,2) = sum(clustmap(:) == cl);
            clustinfo(cl,3) = sum(pos_clustmap(clustmap(:) == cl));
        end
        max_clust_info_pos(p,:) = [clustinfo(find(clustinfo(:,2) == max(clustinfo(:,2)),1,'first'),2) clustinfo(find(clustinfo(:,3) == max(clustinfo(:,3)),1,'first'),3)];
    else
        max_clust_info_pos(p,:) = [0 0];
    end
    clear cl Num clustinfo clustmap
    
    waitbar(p/nPerms,h,sprintf('%d of %d Permutation finished!',p,nPerms))
end
close(h)


%% Z Map thresholded

% Plot significant uncorrected areas
zmapthresh_pos = Real_T;
zmapthresh_pos(zmapthresh_pos < tinv(1-thresh_pval, size(Data,1)-1)) = 0;
[clustmap,Num] = bwlabel(zmapthresh_pos);  %clustinfo = zeros(Num,3);
for cl = 1:Num
    clustinfo(cl,1) = cl;
    clustinfo(cl,2) = sum(clustmap(:) == cl);
    clustinfo(cl,3) = sum(zmapthresh_pos(clustmap(:) == cl));
end
clust_threshold = prctile(max_clust_info_pos(:,2),100-(mcc_cluster_pval)*100);
if(Num > 0)
    for i = 1:size(clustinfo)
        if(clustinfo(i,3) < clust_threshold) 
            zmapthresh_pos(clustmap == clustinfo(i,1)) = 0;
        end
    end
end
clear cl i Num clustmap

zmapthresh_pos(zmapthresh_pos == 0) = NaN;

%% Save into Results

Results.nPerms  = nPerms;
Results.TimeVec = TimeVec;
Results.Hyp_Mat = Hyp_Mat;
Results.Real_T  = Real_T;
Results.Surr_T  = Surr_T;
Results.max_pixel_pvals  = max_pixel_pvals;
Results.max_clust_info_pos  = max_clust_info_pos;
Results.clustinfo = clustinfo;
Results.zmapthresh_pos = zmapthresh_pos;

end