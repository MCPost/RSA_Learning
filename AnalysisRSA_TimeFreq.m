%% RSA with time-frequency data

%% Subject Names
Subj_names = {'AG','CEWD','CM','DS','FVM'};


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


ROI = {'OCC','TMP'};
frqmth = {'FT','WL','FH'};
frqs = {'theta','alpha','beta','gamma1','gamma2'};
blnames = {'NoBL','logBL','pctBL'};
measures = {'LDA','SVM','euclidian'};

for sub = 1:length(Subj_names)
    for fm = frqmth
        curData = load(['RSA_DiffFreq/RSA_Power_',Subj_names{sub}], ['RSA_Power_',fm{1}]);
        for fq = frqs
            for bl = blnames
                for m = 1:length(measures)
                    if(sub == 1)
                        dims = size(curData.(['RSA_Power_',fm{1}]).(ROI{1}).(fq{1}).(bl{1}).RSA_16{1,m});
                        RSA_DiffFreq.(ROI{1}).(fm{1}).(fq{1}).(bl{1}){1,m} = nan(dims(3), dims(1), dims(2));
                        RSA_DiffFreq.(ROI{2}).(fm{1}).(fq{1}).(bl{1}){1,m} = nan(dims(3), dims(1), dims(2));
                    end
                    RSA_DiffFreq.(ROI{1}).(fm{1}).(fq{1}).(bl{1}){1,m}(sub,:,:,:) = permute(curData.(['RSA_Power_',fm{1}]).(ROI{1}).(fq{1}).(bl{1}).RSA_16{1,m},[3 1 2]);
                    RSA_DiffFreq.(ROI{2}).(fm{1}).(fq{1}).(bl{1}){1,m}(sub,:,:,:) = permute(curData.(['RSA_Power_',fm{1}]).(ROI{1}).(fq{1}).(bl{1}).RSA_16{1,m},[3 1 2]);
                end
            end
        end
    end
end
RSA_DiffFreq.TimeVec = curData.(['RSA_Power_',fm{1}]).(ROI{1}).(fq{1}).(bl{1}).TimeVec;

save('RSA_DiffFreq','RSA_DiffFreq','Subj_names','ROI_occ','ROI_temp','ROI','frqmth','frqs','blnames','measures')


