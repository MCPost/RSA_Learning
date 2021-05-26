%% Crosscompare RSA Enc Ret LoadData

function Crosscomparison_RSA_Enc_Ret_LoadData(msr, cfg, pair, save_name)

% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs
ROI = cfg.ROI;

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



%% Load Data 

% Subject Names
load('RSA_Data_Enc','Subj_names')

% Create Data Struct for Encoding
RSA_Data_Enc = [];
RSA_Data_Ret = [];
for sub = 1:length(Subj_names)
    tmp_strct_enc = load(['RSA_Data_wholehead/RSA_Data_Enc_wholehead_',Subj_names{sub}]);
    tmp_strct_ret = load(['RSA_Data_wholehead/RSA_Data_Ret_wholehead_',Subj_names{sub}]);
    if(sub == 1)
        RSA_Data_Enc.Names   = Subj_names;
        for fn = fieldnames(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{1}))'
            RSA_Data_Enc.(fn{1}) = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).(fn{1});
        end
        RSA_Data_Enc = rmfield(RSA_Data_Enc, {'Name','TimeVec1024','RSA_full','MDS_full','RSA_16','Encoding_Data','TrialInfo'...
                                      'rsa_dim','curROI','curROI_name'});
        for r = 1:length(ROI)
            RSA_Data_Enc.([ROI{r},'_ROI']) = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).curROI;
        end

        RSA_Data_Ret.Names   = Subj_names;
        for fn = fieldnames(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{1}))'
            RSA_Data_Ret.(fn{1}) = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).(fn{1});
        end
        RSA_Data_Ret = rmfield(RSA_Data_Ret, {'Name','TimeVec1024','RSA_full','MDS_full','RSA_16','Retrieval_Data','TrialInfo'...
                                      'rsa_dim','curROI','curROI_name'});
        for r = 1:length(ROI)
            RSA_Data_Ret.([ROI{r},'_ROI']) = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).curROI;
        end

    end
    
    for r = 1:length(ROI)
        if(~isempty(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).RSA_full))
            RSA_Data_Enc.(ROI{r}).full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_full{1,msr},[3 1 2]);
        end
        RSA_Data_Enc.(ROI{r}).red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_16{1,msr},[3 1 2]);
        
        if(~isempty(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).RSA_full))
            RSA_Data_Ret.(ROI{r}).full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_full{1,msr},[3 1 2]);
        end
        RSA_Data_Ret.(ROI{r}).red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{r}).RSA_16{1,msr},[3 1 2]);
    end
    
    RSA_Data_Enc.Encoding_Data{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).Encoding_Data;
    RSA_Data_Enc.TrialInfo{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).TrialInfo;

    RSA_Data_Ret.Retrieval_Data{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).Retrieval_Data;
    RSA_Data_Ret.TrialInfo{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).(ROI{1}).TrialInfo;

end
clear tmp_strct_enc tmp_strct_ret



%% Cross-correlate RSA time courses for each subject

%cfg = [];
%cfg.subs             = Subj_names;
%cfg.slide_window     = 0.080;
%cfg.slide_step       = 0.010;
%cfg.fwhm             = 0.040;
%cfg.window_average   = 'gaussian';
%cfg.Hyp_perceptual   = {Perceptual_Mat_full Perceptual_Mat_red16}; %{triu(ones(128)).*~eye(128) triu(ones(16)).*~eye(16)}; %{Perceptual_Mat_full Perceptual_Mat_red16};
%cfg.Hyp_semantic     = {Semantic_Mat_full   Semantic_Mat_red16};
%cfg.ROI              = {'OCC';'TMP';'FRT';'PRT'};
%cfg.timewind1        = [-0.2 1.5];
%cfg.timewind2        = [-2.5 0.2];
%cfg.only16           = true;
%cfg.permtest         = false;
%cfg.n_perms          = 1000;
%cfg.thresh_pval      = 0.05;
%cfg.mcc_cluster_pval = 0.05;
%cfg.ts_os_fac        = 1;
%cfg.matshuffle       = true;
%cfg.studentized      = false;
%Data1 = RSA_Data_Enc; Data2 = RSA_Data_Ret;

if(strcmp(pair,'Enc - Enc'))
    Data1 = RSA_Data_Enc;
    Data2 = RSA_Data_Enc;
elseif(strcmp(pair,'Ret - Ret'))
    Data1 = RSA_Data_Ret;
    Data2 = RSA_Data_Ret;
else
    Data1 = RSA_Data_Enc; 
    Data2 = RSA_Data_Ret;
end

save([save_name,num2str(msr)],'Subj_names')

CrossComp_RSA.(['CrossComp_RSA_',strrep(strrep(measures{msr},' ','_'),'.','')]) = create_xcomp_rsa(cfg, Data1, Data2);

save([save_name,num2str(msr)],'-struct','CrossComp_RSA','-append')

end

