%% RSA Encoding - All Subjects Analysis

%% Subject Names
[~,message,~] = fileattrib('Preproc_EEG_Data\Encoding_object_locked\*');
currentdir = pwd;
filenames = strrep({message([message.directory] == 0).Name}',[currentdir,'\Preproc_EEG_Data\Encoding_object_locked\'],'');
Subj_names = cellfun(@(x) x{1}(10:end), regexp(filenames,'_(\w*).','tokens','once'),'UniformOutput', 0);


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
elecs = {'A1';'A2';'A3';'A4';'A5';'A6';'A7';'A8';'A9';'A10';'A11';'A12';'A13';'A14';'A15';'A16';'A17';...
          'A18';'A19';'A20';'A21';'A22';'A23';'A24';'A25';'A26';'A27';'A28';'A29';'A30';'A31';'A32';'B1';...
          'B2';'B3';'B4';'B5';'B6';'B7';'B8';'B9';'B10';'B11';'B12';'B13';'B14';'B15';'B16';'B17';'B18';...
          'B19';'B20';'B21';'B22';'B23';'B24';'B25';'B26';'B27';'B28';'B29';'B30';'B31';'B32';'C1';'C2';...
          'C3';'C4';'C5';'C6';'C7';'C8';'C9';'C10';'C11';'C12';'C13';'C14';'C15';'C16';'C17';'C18';'C19';...
          'C20';'C21';'C22';'C23';'C24';'C25';'C26';'C27';'C28';'C29';'C30';'C31';'C32';'D1';'D2';'D3';...
          'D4';'D5';'D6';'D7';'D8';'D9';'D10';'D11';'D12';'D13';'D14';'D15';'D16';'D17';'D18';'D19';'D20';...
          'D21';'D22';'D23';'D24';'D25';'D26';'D27';'D28';'D29';'D30';'D31';'D32'};

% All Electrodes
ROI_all_idx = 1:128;

% Occipital Cortex ROI
ROI_occ = {'A8','A9','A10','A11','A12','A13','A14','A15','A16','A17','A21','A22','A23','A24','A25','A26',...
           'A27','A28','A29','A30','B5','B6','B7','B8','B9',...
           'D32','D31','D30','D29','A6','A7','A5','A18','A19','A20','A32','A31','B3','B4','B13','B12','B11','B10'};
ROI_occ_idx = find(cell2mat(cellfun(@(x) any(strcmp(x, ROI_occ)), elecs, 'UniformOutput', 0)));

% Temporal Cortex ROI
ROI_temp = {'B18','B17','B16','B15','B14','B22','B23','B24','B25','B26','B31','B30','B29','B28','B27','C3',...
            'C4','C5','C6','C7','D28','D27','D26','D25','D24','D19','D20','D21','D22','D23','D8','D9','D10',...
            'D11','D12','D7','D6','D5','D4','D3'};
ROI_temp_idx = find(cell2mat(cellfun(@(x) any(strcmp(x, ROI_temp)), elecs, 'UniformOutput', 0)));

%% Measure
msr = 1;


%% Create Data
tmp_strct = load('RSA_Data_Enc');
% RSA_Data_OCC_full.TimeVec = RSA_Data.('RSA_Data_OCC_AG').TimeVec;
% RSA_Data_OCC_16.TimeVec = RSA_Data.('RSA_Data_OCC_AG').TimeVec;
% RSA_Data_TMP_16.TimeVec = RSA_Data.('RSA_Data_TMP_AG').TimeVec;
% if(~isempty(RSA_Data.('RSA_Data_OCC_AG').RSA_full)) 
%     RSA_Data_OCC_full.Data = zeros(size(RSA_Data.('RSA_Data_OCC_AG').RSA_full{1})); RSA_Data_TMP_full.Data = zeros(size(RSA_Data.('RSA_Data_TMP_AG').RSA_full{1})); 
% end
% RSA_Data_OCC_16.Data = zeros(size(RSA_Data.('RSA_Data_OCC_AG').RSA_16{1})); RSA_Data_TMP_16.Data = zeros(size(RSA_Data.('RSA_Data_TMP_AG').RSA_16{1}));

RSA_Data = [];
for sub = 1:length(Subj_names)
    if(sub == 1)
        RSA_Data.Names   = Subj_names;
        for fn = fieldnames(tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC)'
            RSA_Data.(fn{1}) = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.(fn{1});
        end
        RSA_Data = rmfield(RSA_Data, {'Name','TimeVec1024','RSA_full','MDS_full','RSA_16','Encoding_Data','TrialInfo'...
                                      'MDS_16','rsa_dim','mds_dim','mds_error_16','curROI','curROI_name'});
        RSA_Data.Encoding_Data{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.Encoding_Data;
        RSA_Data.TrialInfo{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;
        RSA_Data.OCC_ROI = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.curROI; 
        RSA_Data.TMP_ROI = tmp_strct.(['RSA_Data_',Subj_names{sub}]).TMP.curROI;
    else
        RSA_Data.Encoding_Data{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.Encoding_Data;
        RSA_Data.TrialInfo{sub} = tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;
    end
    
    if(~isempty(tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full))
        RSA_Data.OCC.full_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full{1,msr},[3 1 2]);
        RSA_Data.TMP.full_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_full{1,msr},[3 1 2]);
    end
    RSA_Data.OCC.red16_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_16{1,msr},[3 1 2]);
    RSA_Data.TMP.red16_Data(sub,:,:,:) = permute(tmp_strct.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_16{1,msr},[3 1 2]);
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
%Results = rsa_perm(cfg, RSA_Data.TMP.red16_Data);


% Plot RSA Time series

ROI = {'OCC','TMP'}; r = 1;
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

r1 = 2; r2 = 2;
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
sign_mcc_clust_1 = Results1.zmapthresh_pos;
sign_mcc_clust_1(sign_mcc_clust_1 > 0) = min(get(gca,'ylim'))*0.8;
plot(TimeVec,sign_mcc_clust_1,'bo','MarkerFaceColor','b')
sign_mcc_clust_2 = Results2.zmapthresh_pos;
sign_mcc_clust_2(sign_mcc_clust_2 > 0) = min(get(gca,'ylim'))*0.9;
plot(TimeVec,sign_mcc_clust_2,'ro','MarkerFaceColor','r')
hold off
saveas(gcf,'Results/Enc_TMP_PerceptualvsSemantic_BT-WI_Dat16_LDA.png')
close(gcf)




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