%% Deprecated createRSA function 

%function RSA_Data = createRSA(cfg)

% Parameter Settings
slide_window_s = cfg.slide_window; %s
slide_step_s = cfg.slide_step; %s
if(strcmp(cfg.window_average,'gaussian'))
    average_kern = @(x, dim, gsize) sum(bsxfun(@times, x, (gausswin(gsize)./sum(gausswin(gsize)))'),dim);
    av = 'gaussian'; % Gaussian Average
else
    average_kern = @(x, dim, gsize) mean(x,dim);
    av = 'uniform';  % Uniform Average
end
Name = cfg.Name;
TrialInfo = cfg.TrialInfo;
TimeVec1024 = cfg.TimeVec1024; 
curROI = cfg.curROI;
meas128 = cfg.meas128;
meas16  = cfg.meas16;
MNN = cfg.MNN; 
Cktl_blank_rm = cfg.Cktl_blank_rm;
only16 = cfg.only16;
Data_EEG = cfg.Data;
cfg = rmfield(cfg,'Data');
cfg = rmfield(cfg,'dim');

% Visualization of the Gaussian and Normal (uniform) average weights
%n = 21;
%figure;plot(1:n, gausswin(n)'./sum(gausswin(n)));hold on; plot(1:n, (1/n)*ones(1,n));ylim([0 1.2/sum(gausswin(n))]);xlim([1 n]); hold off


%% Create RSA Matrix
samplingrate = 1024; %Hz
slide_window_idx = round(samplingrate*(slide_window_s/2));
slide_step_idx = round(samplingrate*slide_step_s);
TimeVec_idx = 1:slide_step_idx:length(TimeVec1024);
TimeVec = TimeVec1024(TimeVec_idx);
cfg.TimeVec = TimeVec;


% Save Data
RSA_Data = cfg;
RSA_Data.RSA_full   = []; 
RSA_Data.MDS_full   = []; 

if(~only16)

    fprintf('\n')
    fprintf('Create RSA Matrix (Slide Window: %1.3f s | Slide Step: %1.3f s) with %s window averaging for Subject %s!', slide_window_s, slide_step_s, av, Name)
    fprintf('\n')

    if(strcmp(meas128,'pearson'))
        measures = {'pearson'};
    elseif(strcmp(meas128,'spearman'))
        measures = {'spearman'};
    elseif(strcmp(meas128,'euclidian'))
        measures = {'euclidian'};
        MDS_Mat = {zeros(size(Data_EEG,1),2)};
        MDS_Mat(2,:) = {'euclidian'};
        mds_error_full = {};
        mds_back = zeros(1,length(TimeVec));
    elseif(strcmp(meas128,'mahalanobis'))
        measures = {'mahalanobis'};
        MDS_Mat = {zeros(size(Data_EEG,1),2)};
        MDS_Mat(2,:) = {'mahalanobis'};
        mds_error_full = {};
        mds_back = zeros(1,length(TimeVec));
    elseif(strcmp(meas128,'cosine'))
        measures = {'cosine'};
    else
        measures = {'pearson','spearman','euclidian','mahalanobis','cosine'};
        MDS_Mat = {zeros(size(Data_EEG,1),2), zeros(size(Data_EEG,1),2)};
        MDS_Mat(2,:) = {'euclidian','mahalanobis'};
        mds_error_full = {{},{}};
        mds_back = zeros(2,length(TimeVec));
    end
    RSA_Mat = repmat({repmat(eye(size(Data_EEG,1)),[1 1 length(TimeVec)])},1,length(measures));
    RSA_Mat(2,:) = measures;
    fprintf('\n')
    fprintf('RSA Matrix Computation Progress:  0.0 %%')
    txt = sprintf('RSA Matrix Computation Progress:  0.0 %%');
    for tp = 1:length(TimeVec) 
        time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, size(Data_EEG,3));

        all_patterns = zeros(size(Data_EEG,1), length(curROI));
        for i = 1:size(Data_EEG,1)
            all_patterns(i,:) = average_kern(squeeze(Data_EEG(i,curROI,time_window)), 2, length(time_window))';
        end
        if(Cktl_blank_rm)
            mean_pattern = nanmean(all_patterns,1);
        else
            mean_pattern = zeros(1,size(all_patterns,2));
        end

        for i = 1:size(Data_EEG,1)-1
            for j = (i+1):size(Data_EEG,1)
                cur_trial1 = average_kern(squeeze(Data_EEG(i,curROI,time_window)), 2, length(time_window)) - mean_pattern';
                cur_trial2 = average_kern(squeeze(Data_EEG(j,curROI,time_window)), 2, length(time_window)) - mean_pattern';
                for m = 1:length(measures)
                    RSA_Mat{1,m}(i,j,tp) = dist_meas(cur_trial1, cur_trial2, measures{m});
                end
            end
        end
        
        if(exist('MDS_Mat','var'))
            for m = 1:length(MDS_Mat(2,:))
                if(isempty(find(mds_back(m,:) == 1,1,'last')))
                    try
                        MDS_Mat{1,m}(:,:,tp) = mdscale(RSA_Mat{1,strcmp(RSA_Mat(2,:),MDS_Mat(2,m))}(:,:,tp).*~eye(128) + transpose(RSA_Mat{1,strcmp(RSA_Mat(2,:),MDS_Mat(2,m))}(:,:,tp).*~eye(128)),2,'Start','random');
                        mds_back(m,tp) = 1;
                    catch ME
                        MDS_Mat{1,m}(:,:,tp) = nan(size(MDS_Mat{1,m}(:,:,tp)));
                        mds_error_full{1,m} = [mds_error_full{1,m}; {ME tp}];
                        mds_back(m,tp) = 0;
                    end
                else
                    try
                        MDS_Mat{1,m}(:,:,tp) = mdscale(RSA_Mat{1,strcmp(RSA_Mat(2,:),MDS_Mat(2,m))}(:,:,tp).*~eye(128) + transpose(RSA_Mat{1,strcmp(RSA_Mat(2,:),MDS_Mat(2,m))}(:,:,tp).*~eye(128)),2,'Start',MDS_Mat{1,m}(:,:,find(mds_back(m,:) == 1,1,'last')));
                        mds_back(m,tp) = 1;
                    catch ME
                        MDS_Mat{1,m}(:,:,tp) = nan(size(MDS_Mat{1,m}(:,:,tp)));
                        mds_error_full{1,m} = [mds_error_full{1,m}; {ME tp}];
                        mds_back(m,tp) = 0;
                    end
                end
            end

            %[U,S,V] = svd(MDS_Mat_av(:,:,tp-1)'*MDS_Mat_av(:,:,tp));
            %r = sign(det(V*U'));
            %R = V*[1 0; 0 r]*U';
            %MDS_Mat_av(:,:,tp) = MDS_Mat_av(:,:,tp)*R;
        end

        fprintf(repmat('\b',1,length(txt)))
        fprintf('RSA Matrix Computation Progress: %3.2f %%',(tp/length(TimeVec))*100)
        txt = sprintf('RSA Matrix Computation Progress: %3.2f %%',(tp/length(TimeVec))*100);
    end
    fprintf('\n\n')
    fprintf('Full RSA Matrix for Subject %s created!', Name)
    fprintf('\n')

    % Save Data
    RSA_Data.RSA_full       = RSA_Mat;
    RSA_Data.MDS_full       = MDS_Mat;
    RSA_Data.mds_error_full = mds_error_full;
end
    
% Multivariate Noise Normalization
if(MNN)
    nonnan_ind = find(cell2mat(TrialInfo(:,1)) ~= 0);
    Sigma = zeros(16,length(curROI),length(curROI));
    for c = 1:16
        tmp = zeros(size(Data_EEG,3),length(curROI),length(curROI));
        for t = 1:size(Data_EEG,3)
            dat = Data_EEG(intersect(nonnan_ind, (c-1)*8+(1:8)),curROI,t);
            tmp(t,:,:) = cov1para(dat(~sum(isnan(dat),2),:));
        end
        Sigma(c,:,:) = nanmean(tmp,1);
    end
    tmp = squeeze(mean(Sigma,1));
    sigma_inv = tmp^-0.5;
    Data_EEG_corr = zeros(size(Data_EEG,1), length(curROI), size(Data_EEG,3));
    for t = 1:size(Data_EEG,3)
        Data_EEG_corr(:,:,t) = Data_EEG(:,curROI,t)*sigma_inv;
    end
    fprintf('\n\n')
    fprintf('Multivariate Noise Normalization for Subject %s created!', Name)
    fprintf('\n')
else
    Data_EEG_corr = Data_EEG(:,curROI,:);
end


% Average over Examples: 16 x 16 Matrix
measures = cell(1,length(meas16));
for m = 1:length(meas16)
    if(strcmp(meas16(m),'euclidian'))
        measures(1,m) = {'euclidian c.v.'};
    elseif(strcmp(meas16(m),'euclidian w.c.c.'))
        measures(1,m) = {'euclidian w.c.c.'};
    elseif(strcmp(meas16(m),'pearson'))
        measures(1,m) = {'pearson c.v.'};
    elseif(strcmp(meas16(m),'LDA'))
        measures(1,m) = {'LDA'};
    elseif(strcmp(meas16(m),'SVM'))
        measures(1,m) = {'SVM'};
    elseif(strcmp(meas16(m),'DistCorr'))
        measures(1,m) = {'DistCorr'};
    else
        measures = {'euclidian c.v.','pearson c.v.','LDA','SVM','DistCorr'};
    end
end
RSA_Mat_16x16 = repmat({repmat(zeros(16),[1 1 length(TimeVec)])},1,length(measures));
RSA_Mat_16x16(2,:) = measures;    
MDS_Mat_16x16 = repmat({zeros(16,2,length(TimeVec))},1,length(measures));
MDS_Mat_16x16(2,:) = measures;
mds_error_16 = repmat({{}},1,length(measures));
Data = Data_EEG_corr;
cur_trial1 = zeros(8,length(curROI));
cur_trial2 = zeros(8,length(curROI));
mds_back = zeros(length(measures),length(TimeVec));
fprintf('\n')
fprintf('Progress 16x16:  0.0 %%')
txt = sprintf('Progress 16x16:  0.0 %%');
for tp = 1:length(TimeVec) 
    time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, size(Data,3));
    
    all_patterns = zeros(size(Data,1), length(curROI));
    for i = 1:size(Data,1)
        all_patterns(i,:) = average_kern(squeeze(Data(i,:,time_window)), 2, length(time_window))';
    end
    if(Cktl_blank_rm)
        mean_pattern = nanmean(all_patterns,1);
    else
        mean_pattern = zeros(1,size(all_patterns,2));
    end
    
    ct_x = 1; ct_y = 1; within_dist = zeros(16,1);
    for i = 1:8:size(Data,1)-1
        for j = i:8:size(Data,1)
            if(i ~= j)
                for h = 0:7
                    cur_trial1(h+1,:) = average_kern(squeeze(Data(i+h,:,time_window)), 2, length(time_window))' - mean_pattern;
                    cur_trial2(h+1,:) = average_kern(squeeze(Data(j+h,:,time_window)), 2, length(time_window))' - mean_pattern;
                end
                x = cur_trial1; y = cur_trial2;
                cur_trial1(sum(isnan(cur_trial1),2) > 0,:) = []; cur_trial2(sum(isnan(cur_trial2),2) > 0,:) = []; Cur_trial = [cur_trial1; cur_trial2];
                
                % Create all differences
                poss_folds = [kron(1:size(cur_trial1,1),ones(1,size(cur_trial2,1)))' kron(ones(1,size(cur_trial1,1)),1:size(cur_trial2,1))'];
                diff_crossv = zeros(size(cur_trial1,2), size(poss_folds,1));
                for df = 1:size(poss_folds,1)
                    diff_crossv(:,df) = cur_trial1(poss_folds(df,1),:)' - cur_trial2(poss_folds(df,2),:)';
                end

                % Crossvalidation Folds
                if(sum(strcmp(measures,'euclidian c.v.')) > 0)
                    dist_crossv = zeros(size(diff_crossv,1),2);
                    for fold = 1:size(diff_crossv,2)
                        dist_crossv(fold,1) = diff_crossv(:,fold)'*mean(diff_crossv(:,[1:fold-1 fold+1:end]),2);
                    end
                    RSA_Mat_16x16{1,strcmp(RSA_Mat_16x16(2,:),'euclidian c.v.')}(ct_x,ct_y,tp) = mean(dist_crossv(:,1));
                end
                
                if(sum(strcmp(measures,'euclidian w.c.c.')) > 0)
                    dist_crossv = zeros(size(diff_crossv,1),1);
                    for fold = 1:size(diff_crossv,2)
                        dist_crossv(fold) = diff_crossv(:,fold)'*diff_crossv(:,fold);
                    end
                    if(within_dist(ct_x) == 0)
                        for f1 = 1:size(cur_trial1,1)-1
                            for f2 = f1:size(cur_trial1,1)
                                within_dist(ct_x) = within_dist(ct_x) + (cur_trial1(1,:)-cur_trial1(2,:))*(cur_trial1(1,:)-cur_trial1(2,:))';
                            end
                        end
                    end
                    if(within_dist(ct_y) == 0)
                        for f1 = 1:size(cur_trial2,1)-1
                            for f2 = f1:size(cur_trial2,1)
                                within_dist(ct_y) = within_dist(ct_y) + (cur_trial2(1,:)-cur_trial2(2,:))*(cur_trial2(1,:)-cur_trial2(2,:))';
                            end
                        end
                    end
                    RSA_Mat_16x16{1,strcmp(RSA_Mat_16x16(2,:),'euclidian w.c.c.')}(ct_x,ct_y,tp) = mean(dist_crossv) - sum(within_dist([ct_x ct_y]))/(2*size(cur_trial1,1)*(size(cur_trial2,1)));
                end
                
                if(sum(strcmp(measures,'pearson c.v.')) > 0)
                    dist_crossv = zeros(size(diff_crossv,1),2);
                    for fold = 1:size(diff_crossv,2)
                        dist_crossv(fold,1) = diff_crossv(:,fold)'*mean(diff_crossv(:,[1:fold-1 fold+1:end]),2);
                        %dist_crossv(fold,2) = diff_crossv(:,fold)'*sigma_inv*mean(diff_crossv(:,[1:fold-1 fold+1:end]),2);

                        eps1 = .1 * var(mean(cur_trial1([1:poss_folds(fold,1)-1 poss_folds(fold,1)+1:end],:),1));
                        eps2 = .1 * var(mean(cur_trial2([1:poss_folds(fold,2)-1 poss_folds(fold,2)+1:end],:),1));
                        eps3 = .25 * sqrt(var(mean(cur_trial1([1:poss_folds(fold,1)-1 poss_folds(fold,1)+1:end],:),1)) * var(mean(cur_trial2([1:poss_folds(fold,2)-1 poss_folds(fold,2)+1:end],:),1)));
                        a1 = getfield(cov(cur_trial1(poss_folds(fold,1),:)', mean(cur_trial2([1:poss_folds(fold,2)-1 poss_folds(fold,2)+1:end],:),1)'), {2});
                        a2 = getfield(cov(cur_trial2(poss_folds(fold,2),:)', mean(cur_trial1([1:poss_folds(fold,1)-1 poss_folds(fold,1)+1:end],:),1)'), {2});
                        a3 = max(eps1, getfield(cov(cur_trial1(poss_folds(fold,1),:)', mean(cur_trial1([1:poss_folds(fold,1)-1 poss_folds(fold,1)+1:end],:),1)'), {2}));
                        a4 = max(eps2, getfield(cov(cur_trial2(poss_folds(fold,2),:)', mean(cur_trial2([1:poss_folds(fold,2)-1 poss_folds(fold,2)+1:end],:),1)'), {2}));
                        dist_crossv(fold,2) = 1 - min(1, max(-1, (0.5*(a1 + a2))/max(eps3, sqrt(a3*a4))));
                    end
                    RSA_Mat_16x16{1,strcmp(RSA_Mat_16x16(2,:),'pearson c.v.')}(ct_x,ct_y,tp) = mean(dist_crossv(:,2));
                end
                
                % SVM and LDA
                if(sum(strcmp(measures,'LDA')) > 0)
                    cfg = [];
                    cfg.metric      = 'acc';
                    cfg.cv          = 'holdout';
                    cfg.p           = 0.3;
                    cfg.repeat      = 10;
                    cfg.classifier  = 'lda';
                    cfg.feedback    = 0;
                    RSA_Mat_16x16{1,strcmp(RSA_Mat_16x16(2,:),'LDA')}(ct_x,ct_y,tp) = mv_classify(cfg, Cur_trial, [ones(1,size(cur_trial1,1)) 2*ones(1,size(cur_trial2,1))]');
                end
                if(sum(strcmp(measures,'SVM')) > 0)
                    cfg = [];
                    cfg.metric      = 'acc';
                    cfg.cv          = 'holdout';
                    cfg.p           = 0.3;
                    cfg.repeat      = 10;
                    cfg.classifier  = 'svm';
                    cfg.feedback    = 0;
                    RSA_Mat_16x16{1,strcmp(RSA_Mat_16x16(2,:),'SVM')}(ct_x,ct_y,tp) = mv_classify(cfg, Cur_trial, [ones(1,size(cur_trial1,1)) 2*ones(1,size(cur_trial2,1))]');
                end
                
                % Distance Correlation (Brownian Correlation)
                if(sum(strcmp(meas16,'DistCorr')) > 0)
                    RSA_Mat_16x16{1,strcmp(RSA_Mat_16x16(2,:),'DistCorr')}(ct_x,ct_y,tp) = distCorr(x, y);
                end
            end
            ct_y = ct_y + 1;
        end
        ct_x = ct_x + 1;
        ct_y = ct_x;
    end
    for m = 1:length(measures)
        RSA_Mat_16x16{1,m}(RSA_Mat_16x16{1,m} < 0) = NaN;
        
        if(isempty(find(mds_back(m,:) == 1,1,'last')))
            try
                MDS_Mat_16x16{1,m}(:,:,tp) = mdscale(RSA_Mat_16x16{1,m}(:,:,tp) + RSA_Mat_16x16{1,m}(:,:,tp)',2,'Start','random');
                mds_back(m,tp) = 1;
            catch ME
                MDS_Mat_16x16{1,m}(:,:,tp) = nan(size(MDS_Mat_16x16{1,m}(:,:,tp)));
                mds_error_16{1,m} = [mds_error_16{1,m}; {ME tp}];
                mds_back(m,tp) = 0;
            end
        else
            try
                MDS_Mat_16x16{1,m}(:,:,tp) = mdscale(RSA_Mat_16x16{1,m}(:,:,tp) + RSA_Mat_16x16{1,m}(:,:,tp)',2,'Start',MDS_Mat_16x16{1,m}(:,:,find(mds_back(m,:) == 1,1,'last')));
                mds_back(m,tp) = 1;
            catch ME
                MDS_Mat_16x16{1,m}(:,:,tp) = nan(size(MDS_Mat_16x16{1,m}(:,:,tp)));
                mds_error_16{1,m} = [mds_error_16{1,m}; {ME tp}];
                mds_back(m,tp) = 0;
            end
        end
       
        %[U,S,V] = svd(MDS_Mat_av(:,:,tp-1)'*MDS_Mat_av(:,:,tp));
        %r = sign(det(V*U'));
        %R = V*[1 0; 0 r]*U';
        %MDS_Mat_av(:,:,tp) = MDS_Mat_av(:,:,tp)*R;
        
    end
    
    
    fprintf(repmat('\b',1,length(txt)))
    fprintf('Progress 16x16: %3.2f %%',(tp/length(TimeVec))*100)
    txt = sprintf('Progress 16x16: %3.2f %%',(tp/length(TimeVec))*100);
end
fprintf('\n\n')
fprintf('16x16 RSA Matrix for Subject %s created!', Name)
fprintf('\n')


%% Save Data
RSA_Data.RSA_16         = RSA_Mat_16x16; 
RSA_Data.MDS_16         = MDS_Mat_16x16;
RSA_Data.rsa_dim        = 'trl_trl_time';
RSA_Data.mds_dim        = 'trl_component_time';
RSA_Data.mds_error_16   = mds_error_16;


%end