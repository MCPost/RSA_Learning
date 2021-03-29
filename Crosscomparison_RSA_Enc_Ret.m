%% Cross-comparison of RSA time courses

% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs

% Measure
msr = 1;


%% Load Data 

% Encoding

% Subject Names
load('RSA_Data_Enc','Subj_names')

% Create Data Struct for Encoding
tmp_strct_enc = load('RSA_Data_Enc');
tmp_strct_ret = load('RSA_Data_Ret');
RSA_Data_Enc = [];
RSA_Data_Ret = [];
for sub = 1:length(Subj_names)
    if(sub == 1)
        RSA_Data_Enc.Names   = Subj_names;
        for fn = fieldnames(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC)'
            RSA_Data_Enc.(fn{1}) = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.(fn{1});
        end
        RSA_Data_Enc = rmfield(RSA_Data_Enc, {'Name','TimeVec1024','RSA_full','MDS_full','RSA_16','Encoding_Data','TrialInfo'...
                                      'rsa_dim','curROI','curROI_name'});
        RSA_Data_Enc.Encoding_Data{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.Encoding_Data;
        RSA_Data_Enc.TrialInfo{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;
        RSA_Data_Enc.OCC_ROI = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.curROI; 
        RSA_Data_Enc.TMP_ROI = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).TMP.curROI;
        RSA_Data_Enc.FRT_ROI = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).FRT.curROI;
        RSA_Data_Enc.CNT_ROI = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).CNT.curROI;
        RSA_Data_Enc.PRT_ROI = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).PRT.curROI;
        
        RSA_Data_Ret.Names   = Subj_names;
        for fn = fieldnames(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC)'
            RSA_Data_Ret.(fn{1}) = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.(fn{1});
        end
        RSA_Data_Ret = rmfield(RSA_Data_Ret, {'Name','TimeVec1024','RSA_full','MDS_full','RSA_16','Retrieval_Data','TrialInfo'...
                                      'rsa_dim','curROI','curROI_name'});
        RSA_Data_Ret.Retrieval_Data{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.Retrieval_Data;
        RSA_Data_Ret.TrialInfo{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;
        RSA_Data_Ret.OCC_ROI = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.curROI; 
        RSA_Data_Ret.TMP_ROI = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).TMP.curROI;
        RSA_Data_Ret.FRT_ROI = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).FRT.curROI;
        RSA_Data_Ret.CNT_ROI = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).CNT.curROI;
        RSA_Data_Ret.PRT_ROI = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).PRT.curROI;
        
        
    else
        RSA_Data_Enc.Encoding_Data{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.Encoding_Data;
        RSA_Data_Enc.TrialInfo{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;
        
        RSA_Data.Retrieval_Data{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.Retrieval_Data;
        RSA_Data.TrialInfo{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;
    end
    
    if(~isempty(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full))
        RSA_Data_Enc.OCC.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Enc.TMP.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Enc.FRT.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Enc.CNT.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).CNT.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Enc.PRT.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_full{1,msr},[3 1 2]);
    end
    RSA_Data_Enc.OCC.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Enc.TMP.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Enc.FRT.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Enc.CNT.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).CNT.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Enc.PRT.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_16{1,msr},[3 1 2]);
    
    if(~isempty(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full))
        RSA_Data_Ret.OCC.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Ret.TMP.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Ret.FRT.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Ret.CNT.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).CNT.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Ret.PRT.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_full{1,msr},[3 1 2]);
    end
    RSA_Data_Ret.OCC.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Ret.TMP.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Ret.FRT.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Ret.CNT.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).CNT.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Ret.PRT.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_16{1,msr},[3 1 2]);
end
clear tmp_strct_enc tmp_strct_ret

% Hypotheses Matrix
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




%% Cross-correlate RSA time courses for each subject

ROI = {'OCC','TMP','FRT','CNT','PRT'};

sub = 15;

per_ind_wi = Perceptual_Mat_red16(:) > 0;
per_ind_bt = Perceptual_Mat_red16(:) < 0;
sem_ind_wi = Semantic_Mat_red16(:) > 0;
sem_ind_bt = Semantic_Mat_red16(:) < 0;

CrossComp_RSA = [];
CrossComp_RSA.OCC = zeros(length(Subj_names),4,length(RSA_Data_Ret.TimeVec),length(RSA_Data_Enc.TimeVec));
CrossComp_RSA.TMP = zeros(length(Subj_names),4,length(RSA_Data_Ret.TimeVec),length(RSA_Data_Enc.TimeVec));
CrossComp_RSA.FRT = zeros(length(Subj_names),4,length(RSA_Data_Ret.TimeVec),length(RSA_Data_Enc.TimeVec));
CrossComp_RSA.CNT = zeros(length(Subj_names),4,length(RSA_Data_Ret.TimeVec),length(RSA_Data_Enc.TimeVec));
CrossComp_RSA.PRT = zeros(length(Subj_names),4,length(RSA_Data_Ret.TimeVec),length(RSA_Data_Enc.TimeVec));
fprintf('\n')
nbytes = fprintf('Progress Cross Correlation:  0.0 %%');
for tp_enc = 1:length(RSA_Data_Enc.TimeVec)
    
    for tp_ret = 1:length(RSA_Data_Ret.TimeVec)
        
        CrossComp_RSA.OCC(:,1,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{1}).red16_Data(:,tp_enc,per_ind_wi))', squeeze(RSA_Data_Ret.(ROI{1}).red16_Data(:,tp_ret,per_ind_wi))', 'type', 'spearman'));
        CrossComp_RSA.OCC(:,2,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{1}).red16_Data(:,tp_enc,per_ind_bt))', squeeze(RSA_Data_Ret.(ROI{1}).red16_Data(:,tp_ret,per_ind_bt))', 'type', 'spearman'));
        CrossComp_RSA.OCC(:,3,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{1}).red16_Data(:,tp_enc,sem_ind_wi))', squeeze(RSA_Data_Ret.(ROI{1}).red16_Data(:,tp_ret,sem_ind_wi))', 'type', 'spearman'));
        CrossComp_RSA.OCC(:,4,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{1}).red16_Data(:,tp_enc,sem_ind_bt))', squeeze(RSA_Data_Ret.(ROI{1}).red16_Data(:,tp_ret,sem_ind_bt))', 'type', 'spearman'));
        
        CrossComp_RSA.TMP(:,1,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{2}).red16_Data(:,tp_enc,per_ind_wi))', squeeze(RSA_Data_Ret.(ROI{2}).red16_Data(:,tp_ret,per_ind_wi))', 'type', 'spearman'));
        CrossComp_RSA.TMP(:,2,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{2}).red16_Data(:,tp_enc,per_ind_bt))', squeeze(RSA_Data_Ret.(ROI{2}).red16_Data(:,tp_ret,per_ind_bt))', 'type', 'spearman'));
        CrossComp_RSA.TMP(:,3,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{2}).red16_Data(:,tp_enc,sem_ind_wi))', squeeze(RSA_Data_Ret.(ROI{2}).red16_Data(:,tp_ret,sem_ind_wi))', 'type', 'spearman'));
        CrossComp_RSA.TMP(:,4,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{2}).red16_Data(:,tp_enc,sem_ind_bt))', squeeze(RSA_Data_Ret.(ROI{2}).red16_Data(:,tp_ret,sem_ind_bt))', 'type', 'spearman'));
        
        CrossComp_RSA.FRT(:,1,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{3}).red16_Data(:,tp_enc,per_ind_wi))', squeeze(RSA_Data_Ret.(ROI{3}).red16_Data(:,tp_ret,per_ind_wi))', 'type', 'spearman'));
        CrossComp_RSA.FRT(:,2,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{3}).red16_Data(:,tp_enc,per_ind_bt))', squeeze(RSA_Data_Ret.(ROI{3}).red16_Data(:,tp_ret,per_ind_bt))', 'type', 'spearman'));
        CrossComp_RSA.FRT(:,3,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{3}).red16_Data(:,tp_enc,sem_ind_wi))', squeeze(RSA_Data_Ret.(ROI{3}).red16_Data(:,tp_ret,sem_ind_wi))', 'type', 'spearman'));
        CrossComp_RSA.FRT(:,4,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{3}).red16_Data(:,tp_enc,sem_ind_bt))', squeeze(RSA_Data_Ret.(ROI{3}).red16_Data(:,tp_ret,sem_ind_bt))', 'type', 'spearman'));
        
        CrossComp_RSA.CNT(:,1,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{4}).red16_Data(:,tp_enc,per_ind_wi))', squeeze(RSA_Data_Ret.(ROI{4}).red16_Data(:,tp_ret,per_ind_wi))', 'type', 'spearman'));
        CrossComp_RSA.CNT(:,2,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{4}).red16_Data(:,tp_enc,per_ind_bt))', squeeze(RSA_Data_Ret.(ROI{4}).red16_Data(:,tp_ret,per_ind_bt))', 'type', 'spearman'));
        CrossComp_RSA.CNT(:,3,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{4}).red16_Data(:,tp_enc,sem_ind_wi))', squeeze(RSA_Data_Ret.(ROI{4}).red16_Data(:,tp_ret,sem_ind_wi))', 'type', 'spearman'));
        CrossComp_RSA.CNT(:,4,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{4}).red16_Data(:,tp_enc,sem_ind_bt))', squeeze(RSA_Data_Ret.(ROI{4}).red16_Data(:,tp_ret,sem_ind_bt))', 'type', 'spearman'));
        
        CrossComp_RSA.PRT(:,1,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{5}).red16_Data(:,tp_enc,per_ind_wi))', squeeze(RSA_Data_Ret.(ROI{5}).red16_Data(:,tp_ret,per_ind_wi))', 'type', 'spearman'));
        CrossComp_RSA.PRT(:,2,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{5}).red16_Data(:,tp_enc,per_ind_bt))', squeeze(RSA_Data_Ret.(ROI{5}).red16_Data(:,tp_ret,per_ind_bt))', 'type', 'spearman'));
        CrossComp_RSA.PRT(:,3,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{5}).red16_Data(:,tp_enc,sem_ind_wi))', squeeze(RSA_Data_Ret.(ROI{5}).red16_Data(:,tp_ret,sem_ind_wi))', 'type', 'spearman'));
        CrossComp_RSA.PRT(:,4,tp_ret,tp_enc) = diag(corr(squeeze(RSA_Data_Enc.(ROI{5}).red16_Data(:,tp_enc,sem_ind_bt))', squeeze(RSA_Data_Ret.(ROI{5}).red16_Data(:,tp_ret,sem_ind_bt))', 'type', 'spearman'));
        
        %CrossComp_RSA_Mat(tp_ret,tp_enc) = fast_corr(squeeze(RSA_Data_Enc.(ROI{1}).red16_Data(sub,tp_enc,mat_ind)), squeeze(RSA_Data_Ret.(ROI{1}).red16_Data(sub,tp_ret,mat_ind)));
        %CrossComp_RSA_Mat(1,tp_ret,tp_enc) = corr(squeeze(RSA_Data_Enc.(ROI{2}).red16_Data(sub,tp_enc,mat_ind_wi)), squeeze(RSA_Data_Ret.(ROI{2}).red16_Data(sub,tp_ret,mat_ind_wi)), 'type', 'spearman');
        %CrossComp_RSA_Mat(2,tp_ret,tp_enc) = corr(squeeze(RSA_Data_Enc.(ROI{2}).red16_Data(sub,tp_enc,mat_ind_bt)), squeeze(RSA_Data_Ret.(ROI{2}).red16_Data(sub,tp_ret,mat_ind_bt)), 'type', 'spearman');
        
        fprintf(repmat('\b',1,nbytes))
        nbytes = fprintf('Progress Cross Correlation: %3.2f %%',((tp_ret + length(RSA_Data_Ret.TimeVec)*(tp_enc - 1)) / (length(RSA_Data_Enc.TimeVec)*length(RSA_Data_Ret.TimeVec)))*100);
    end
    
end
fprintf('\n')
CrossComp_RSA.OCC = flip(CrossComp_RSA.OCC,3);
CrossComp_RSA.TMP = flip(CrossComp_RSA.TMP,3);
CrossComp_RSA.FRT = flip(CrossComp_RSA.FRT,3);
CrossComp_RSA.CNT = flip(CrossComp_RSA.CNT,3);
CrossComp_RSA.PRT = flip(CrossComp_RSA.PRT,3);


load('CrossComp_RSA')


figure
subplot(1,2,1)
contourf(RSA_Data_Enc.TimeVec, RSA_Data_Ret.TimeVec, squeeze(nanmean(CrossComp_RSA.OCC(:,2,:,:) - CrossComp_RSA.OCC(:,1,:,:),1)), 40,'linestyle','none')
caxis([-.06 .06]); set(gca,'xlim', [-0.2 1], 'ylim', [-3 0.2])
hold on
plot([0 0],[RSA_Data_Ret.TimeVec(1) RSA_Data_Ret.TimeVec(end)],'--k','linewidth',1.5)
plot([RSA_Data_Enc.TimeVec(1) RSA_Data_Enc.TimeVec(end)],[0 0],'--k','linewidth',1.5)
hold off
subplot(1,2,2)
contourf(RSA_Data_Enc.TimeVec, RSA_Data_Ret.TimeVec, squeeze(nanmean(CrossComp_RSA.TMP(:,4,:,:) - CrossComp_RSA.TMP(:,3,:,:),1)), 40,'linestyle','none')
caxis([-.06 .06]); set(gca,'xlim', [-0.2 1], 'ylim', [-3 0.2])
hold on
plot([0 0],[RSA_Data_Ret.TimeVec(1) RSA_Data_Ret.TimeVec(end)],'--k','linewidth',1.5)
plot([RSA_Data_Enc.TimeVec(1) RSA_Data_Enc.TimeVec(end)],[0 0],'--k','linewidth',1.5)
hold off




