%% Load Data and Baseline Correction and/or Artifact Rejection

function Data_EEG = preproc_EEGData(cfg, Data_Struct)


Name = cfg.Name;
chan_idx = cfg.chan_idx; 
Art_corr = cfg.Art_corr; 
BL_corr = cfg.BL_corr; 
BL_wind = cfg.BL_wind;
Data = Data_Struct.Data;
cfg.Encoding_Data = Data_Struct.Encoding_Data; 
cfg.TrialInfo = Data_Struct.TrialInfo;
cfg.TimeVec1024 =  Data_Struct.TimeVec1024;

if(Art_corr)
    fprintf('\n')
    fprintf('Artefact Rejection for Subject %s!', Name)
    fprintf('\n')
    load ArtifactRej
    for i = 1:size(trial2rm.(Name),1)
        if(~strcmp(trial2rm.(Name){i,2},'all'))
            Data(trial2rm.(Name){i,2},trial2rm.(Name){i,1},:) = NaN;
        end
    end
    cfg.chan_idx_usable = chan_idx(~ismember(chan_idx, cell2mat(trial2rm.(Name)(strcmp(trial2rm.(Name)(:,2),'all'),1))));
end

% Baseline Correction
if(~strcmp(BL_corr,'no'))
    [~,BL_wind_idx] = min(abs(bsxfun(@minus, [Data_Struct.TimeVec1024; Data_Struct.TimeVec1024], BL_wind')),[],2);
    BL_mean = mean(Data(:,:,BL_wind_idx(1):BL_wind_idx(2)),3);
    if(strcmp(BL_corr,'covBL'))
        for ch = 1:128
            na_idx = find(~sum(isnan(BL_mean(:,ch)),2));
            X = BL_mean(na_idx,ch); %[ones(size(BL_mean,1),1) BL_mean(:,ch)];
            Data(na_idx,ch,:) = squeeze(Data(na_idx,ch,:)) - X*((X'*X)\X'*squeeze(Data(na_idx,ch,:)));
        end
    elseif(strcmp(BL_corr,'minBL'))
        Data = bsxfun(@minus, Data, BL_mean);
    end
end

% Save
Data_EEG = cfg;
Data_EEG.Data = Data;
Data_EEG.dim = 'trl_chan_time';

end
