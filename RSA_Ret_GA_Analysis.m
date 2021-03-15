%% RSA Retrieval - All Subjects Analysis

%% Subject Names
[~,message,~] = fileattrib('Preproc_EEG_Data\Retrieval_response_locked\*');
currentdir = pwd;
filenames = strrep({message([message.directory] == 0).Name}',[currentdir,'\Preproc_EEG_Data\Retrieval_response_locked\'],'');
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
tmp_strct = load('RSA_Data_Ret');
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
