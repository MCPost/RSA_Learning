%% Create Source Level Data

%% Subject Names
[~,message,~] = fileattrib('Preproc_EEG_Data/Encoding_object_locked/*');
currentdir = pwd;
filenames = strrep({message([message.directory] == 0).Name}',[currentdir,'\Preproc_EEG_Data\Encoding_object_locked\'],'');
Subj_names = cellfun(@(x) x{1}(10:end), regexp(filenames,'_(\w*).mat','tokens','once'),'UniformOutput', 0);

%% Clear RW (Sub 20) because not enough trials
Subj_names(strcmp(Subj_names,'RW')) = [];

%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs

load Source_Location_Headmodels/vol
load Source_Location_Headmodels/grid
load Source_Location_Headmodels/mri

for sub = 1:length(Subj_names)
    
    %% Import Data
    load(['Preproc_EEG_Data/Retrieval_response_locked/Ret_EEG_Data_',Subj_names{sub}]);
    eeg_cfg = [];
    eeg_cfg.Name = Subj_names{sub};
    eeg_cfg.chan_label = elecs;
    eeg_cfg.chan_idx = 1:128; 
    eeg_cfg.Art_corr = false;
    eeg_cfg.BL_corr = 'yes';
    eeg_cfg.BL_wind = [-4 0.2];
    [Data_EEG] = preproc_EEGData(eeg_cfg, Ret_Data_EEG);
    
    % Get rid of nan trials
    Data_EEG.Data(cell2mat(Data_EEG.TrialInfo(:,1)) == 0,:,:) = [];
    Data_EEG.TrialInfo(cell2mat(Data_EEG.TrialInfo(:,1)) == 0,:) = [];
    
    % Get into Fieldtrip Format
    % avg, var, time, dof, label, trial, dimord, elec, trialinfo
    dataset = [];
    dataset.label = Data_EEG.chan_label;
    dataset.time = Data_EEG.TimeVec1024;
    dataset.trial = Data_EEG.Data;
    dataset.dimord = 'rpt_chan_time';
    dataset.trialinfo = Data_EEG.TrialInfo;
    
    % timelock data
    cfg                   = [];
    cfg.keeptrials        = 'yes';
    cfg.covariance        = 'yes';
    cfg.covariancewindow  = 'all';
    timelock              = ft_timelockanalysis(cfg, dataset);
    
    load(['Source_Location_Headmodels/ElectrodeMri/elec_aligned_',Subj_names{sub},'.mat']);
    eval(['elec = elec_aligned_',Subj_names{sub}]);
    
    elec.chanpos = elec.chanpos(1:end-3,:);
    elec.elecpos = elec.elecpos(1:end-3,:);
    elec.label = timelock.label;
    
    cfg = [];
    cfg.sourcemodel = grid;     %% where are the sources?
    cfg.headmodel   = vol;      %% how do currents spread?
    cfg.elec        = elec;     %% where are the sensors?

    % how do sources and sensors connect?
    sourcemodel_and_leadfield = ft_prepare_leadfield(cfg);

    
    
    cfg                   = [];
    cfg.method            = 'lcmv';
    cfg.sourcemodel       = sourcemodel_and_leadfield;
    cfg.headmodel         = vol;
    cfg.elec              = elec;
    cfg.senstype          = 'EEG';
    cfg.keeptrials        = 'yes';
    cfg.lcmv.keepfilter   = 'yes';
    cfg.lcmv.fixedori     = 'yes';
    cfg.lcmv.lambda       = 5; % 5% regularisation
    %cfg.normalize         = 'yes';
    source_tmp            = ft_sourceanalysis(cfg, timelock);
    filters               = cell2mat(source_tmp.avg.filter(source_tmp.inside)); % keep filters
    
    % combine filters with raw data to "beam" scalp data to source space
    source            = [];
    %     source.sampleinfo = dataset.sampleinfo;
    source.time       = dataset.time;
    source.trialinfo  = dataset.trialinfo;
    source.fsample    = 256;
    
    % create labels for each virtual electrode
    for c = 1 : sum(grid.inside)
        label{c,1} = ['S' num2str(c)];
    end
    
    source.label = label;
    
    % for each trial, apply filters to the recorded data
    for j = 1:size(dataset.trial,1)
        source.trial{1,j} = single(filters*squeeze(dataset.trial(j,:,:)));
    end
    
end



