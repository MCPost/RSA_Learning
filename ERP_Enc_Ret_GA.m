%% Event-related Potentials Encoding and Retrieval

%% Subject Names
[~,message,~] = fileattrib('Preproc_EEG_Data/Encoding_object_locked/*');
currentdir = pwd;
filenames = strrep({message([message.directory] == 0).Name}',[currentdir,'\Preproc_EEG_Data\Encoding_object_locked\'],'');
Subj_names = cellfun(@(x) x{1}(10:end), regexp(filenames,'_(\w*).mat','tokens','once'),'UniformOutput', 0);


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs


% Load example struct
exmpl_str = load('Prepro_retrieval_RL_4seconds_full_DS', 'r_1_DS');
exmpl_str = exmpl_str.r_1_DS;

%% Encoding

Enc_SubData_Perceptual = cell(2,length(Subj_names));
Enc_SubData_Semantic = cell(2,length(Subj_names));
Enc_SubData_Category = cell(8,length(Subj_names));
for sub = 1:length(Subj_names)
    
    %% Import Data
    load(['Preproc_EEG_Data/Encoding_object_locked/Enc_EEG_Data_',Subj_names{sub}])
    substruct        = [];
    substruct.time   = Enc_Data_EEG.TimeVec1024;
    substruct.label  = exmpl_str.label;
    substruct.avg    = [];
    substruct.dimord = 'chan_time'; %exmpl_str.dimord;
    substruct.elec   = exmpl_str.elec;
    
    % Perceptual
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,5)) == 1,:,:),1));
    Enc_SubData_Perceptual{1,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,5)) == 2,:,:),1));
    Enc_SubData_Perceptual{2,sub} = substruct;
    
    % Semantic
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,7)) == 1,:,:),1));
    Enc_SubData_Semantic{1,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,7)) == 2,:,:),1));
    Enc_SubData_Semantic{2,sub} = substruct;
    
    % Category
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,9)) == 1,:,:),1));
    Enc_SubData_Category{1,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,9)) == 2,:,:),1));
    Enc_SubData_Category{2,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,9)) == 3,:,:),1));
    Enc_SubData_Category{3,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,9)) == 4,:,:),1));
    Enc_SubData_Category{4,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,9)) == 5,:,:),1));
    Enc_SubData_Category{5,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,9)) == 6,:,:),1));
    Enc_SubData_Category{6,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,9)) == 7,:,:),1));
    Enc_SubData_Category{7,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Enc_Data_EEG.Data(cell2mat(Enc_Data_EEG.TrialInfo(:,9)) == 8,:,:),1));
    Enc_SubData_Category{8,sub} = substruct;
    
    clear Enc_Data_EEG
end


% Grand Average ERPs for Encoding

% Perceptual
cfg = [];
[Enc_GA_Perc_Drawing] = ft_timelockgrandaverage(cfg, Enc_SubData_Perceptual{1,:});
[Enc_GA_Perc_Picture] = ft_timelockgrandaverage(cfg, Enc_SubData_Perceptual{2,:});

cfg = [];
cfg.parameter = 'avg';
cfg.baseline = [-0.2 0];
cfg.layout = 'biosemi128.lay';
cfg.showlegend = 'yes';
ft_multiplotER(cfg, Enc_GA_Perc_Drawing, Enc_GA_Perc_Picture)

cfg = [];
cfg.parameter = 'avg';
cfg.baseline = [-0.2 0];
cfg.channel = ROI_occ';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Enc_GA_Perc_Drawing, Enc_GA_Perc_Picture)
cfg.channel = ROI_temp';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Enc_GA_Perc_Drawing, Enc_GA_Perc_Picture)

% Semantic
cfg = [];
[Enc_GA_Sem_Animate] = ft_timelockgrandaverage(cfg, Enc_SubData_Semantic{1,:});
[Enc_GA_Sem_NonAnimate] = ft_timelockgrandaverage(cfg, Enc_SubData_Semantic{2,:});

cfg = [];
cfg.parameter = 'avg';
cfg.baseline = [-0.2 0];
cfg.layout = 'biosemi128.lay';
cfg.showlegend = 'yes';
ft_multiplotER(cfg, Enc_GA_Sem_Animate, Enc_GA_Sem_NonAnimate)

cfg = [];
cfg.parameter = 'avg';
cfg.baseline = [-0.2 0];
cfg.channel = ROI_occ';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Enc_GA_Sem_Animate, Enc_GA_Sem_NonAnimate)
cfg.channel = ROI_temp';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Enc_GA_Sem_Animate, Enc_GA_Sem_NonAnimate)

% Category
cfg = [];
[Enc_GA_Cat_Insect] = ft_timelockgrandaverage(cfg, Enc_SubData_Category{1,:});
[Enc_GA_Cat_Bird] = ft_timelockgrandaverage(cfg, Enc_SubData_Category{2,:});
[Enc_GA_Cat_Mammal] = ft_timelockgrandaverage(cfg, Enc_SubData_Category{3,:});
[Enc_GA_Cat_Sea] = ft_timelockgrandaverage(cfg, Enc_SubData_Category{4,:});
[Enc_GA_Cat_Fruit] = ft_timelockgrandaverage(cfg, Enc_SubData_Category{5,:});
[Enc_GA_Cat_Electronic] = ft_timelockgrandaverage(cfg, Enc_SubData_Category{6,:});
[Enc_GA_Cat_Veggie] = ft_timelockgrandaverage(cfg, Enc_SubData_Category{7,:});
[Enc_GA_Cat_Clothe] = ft_timelockgrandaverage(cfg, Enc_SubData_Category{8,:});


cfg = [];
cfg.parameter = 'avg';
cfg.baseline = [-0.2 0];
cfg.channel = ROI_occ';
cfg.showlegend = 'yes';
cfg.linestyle = {'-','-','-','-','--','--','--','--'};
cfg.linecolor = [0, 0.4470, 0.7410; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.6350, 0.0780, 0.1840;...
                 0, 0.4470, 0.7410; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.6350, 0.0780, 0.1840];
cfg.linewidth = 1.5;
ft_singleplotER(cfg, Enc_GA_Cat_Insect, Enc_GA_Cat_Bird, Enc_GA_Cat_Mammal, Enc_GA_Cat_Sea, Enc_GA_Cat_Fruit, Enc_GA_Cat_Electronic, Enc_GA_Cat_Veggie, Enc_GA_Cat_Clothe)
cfg.channel = ROI_temp';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Enc_GA_Cat_Insect, Enc_GA_Cat_Bird, Enc_GA_Cat_Mammal, Enc_GA_Cat_Sea, Enc_GA_Cat_Fruit, Enc_GA_Cat_Electronic, Enc_GA_Cat_Veggie, Enc_GA_Cat_Clothe)





%% Retrieval

Ret_SubData_Perceptual = cell(2,length(Subj_names));
Ret_SubData_Semantic = cell(2,length(Subj_names));
Ret_SubData_Category = cell(8,length(Subj_names));
for sub = 1:length(Subj_names)
    
    %% Import Data
    load(['Preproc_EEG_Data/Retrieval_response_locked/Ret_EEG_Data_',Subj_names{sub}])
    eeg_cfg = [];
    eeg_cfg.Name = Subj_names{sub};
    eeg_cfg.chan_label = elecs;
    eeg_cfg.chan_idx = ROI_all_idx'; 
    eeg_cfg.Art_corr = false;
    eeg_cfg.BL_corr = 'demean';
    eeg_cfg.BL_wind = [-0.26 -0.05];
    [Data_EEG] = preproc_EEGData(eeg_cfg, Ret_Data_EEG);
    Ret_Data_EEG.Data = Data_EEG.Data;
    substruct        = [];
    substruct.time   = Ret_Data_EEG.TimeVec1024;
    substruct.label  = exmpl_str.label;
    substruct.avg    = [];
    substruct.dimord = 'chan_time'; %exmpl_str.dimord;
    substruct.elec   = exmpl_str.elec;
    
    % Perceptual
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,5)) == 1,:,:),1));
    Ret_SubData_Perceptual{1,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,5)) == 2,:,:),1));
    Ret_SubData_Perceptual{2,sub} = substruct;
    
    % Semantic
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,7)) == 1,:,:),1));
    Ret_SubData_Semantic{1,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,7)) == 2,:,:),1));
    Ret_SubData_Semantic{2,sub} = substruct;
    
    % Category
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,9)) == 1,:,:),1));
    Ret_SubData_Category{1,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,9)) == 2,:,:),1));
    Ret_SubData_Category{2,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,9)) == 3,:,:),1));
    Ret_SubData_Category{3,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,9)) == 4,:,:),1));
    Ret_SubData_Category{4,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,9)) == 5,:,:),1));
    Ret_SubData_Category{5,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,9)) == 6,:,:),1));
    Ret_SubData_Category{6,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,9)) == 7,:,:),1));
    Ret_SubData_Category{7,sub} = substruct;
    substruct.avg    = squeeze(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,9)) == 8,:,:),1));
    Ret_SubData_Category{8,sub} = substruct;
    
    clear Ret_Data_EEG Data_EEG
end


% Grand Average ERPs for Retrieval

% Perceptual
cfg = [];
[Ret_GA_Perc_Drawing] = ft_timelockgrandaverage(cfg, Ret_SubData_Perceptual{1,:});
[Ret_GA_Perc_Picture] = ft_timelockgrandaverage(cfg, Ret_SubData_Perceptual{2,:});

cfg = [];
cfg.parameter = 'avg';
%cfg.baseline = [-0.2 0];
cfg.layout = 'biosemi128.lay';
cfg.showlegend = 'yes';
ft_multiplotER(cfg, Ret_GA_Perc_Drawing, Ret_GA_Perc_Picture)

cfg = [];
cfg.parameter = 'avg';
%cfg.baseline = [-0.2 0];
cfg.channel = ROI_occ';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Ret_GA_Perc_Drawing, Ret_GA_Perc_Picture)
cfg.channel = ROI_temp';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Ret_GA_Perc_Drawing, Ret_GA_Perc_Picture)

% Semantic
cfg = [];
[Ret_GA_Sem_Animate] = ft_timelockgrandaverage(cfg, Ret_SubData_Semantic{1,:});
[Ret_GA_Sem_NonAnimate] = ft_timelockgrandaverage(cfg, Ret_SubData_Semantic{2,:});

cfg = [];
cfg.parameter = 'avg';
%cfg.baseline = [-0.2 0];
cfg.ylim = [-2 2];
cfg.layout = 'biosemi128.lay';
cfg.showlegend = 'yes';
ft_multiplotER(cfg, Ret_GA_Sem_Animate, Ret_GA_Sem_NonAnimate)

cfg = [];
cfg.parameter = 'avg';
%cfg.baseline = [-0.2 0];
cfg.channel = ROI_occ';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Ret_GA_Sem_Animate, Ret_GA_Sem_NonAnimate)
cfg.channel = ROI_temp';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Ret_GA_Sem_Animate, Ret_GA_Sem_NonAnimate)

% Category
cfg = [];
[Ret_GA_Cat_Insect] = ft_timelockgrandaverage(cfg, Ret_SubData_Category{1,:});
[Ret_GA_Cat_Bird] = ft_timelockgrandaverage(cfg, Ret_SubData_Category{2,:});
[Ret_GA_Cat_Mammal] = ft_timelockgrandaverage(cfg, Ret_SubData_Category{3,:});
[Ret_GA_Cat_Sea] = ft_timelockgrandaverage(cfg, Ret_SubData_Category{4,:});
[Ret_GA_Cat_Fruit] = ft_timelockgrandaverage(cfg, Ret_SubData_Category{5,:});
[Ret_GA_Cat_Electronic] = ft_timelockgrandaverage(cfg, Ret_SubData_Category{6,:});
[Ret_GA_Cat_Veggie] = ft_timelockgrandaverage(cfg, Ret_SubData_Category{7,:});
[Ret_GA_Cat_Clothe] = ft_timelockgrandaverage(cfg, Ret_SubData_Category{8,:});


cfg = [];
cfg.parameter = 'avg';
%cfg.baseline = [-0.2 0];
cfg.channel = ROI_occ';
cfg.showlegend = 'yes';
cfg.linestyle = {'-','-','-','-','--','--','--','--'};
cfg.linecolor = [0, 0.4470, 0.7410; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.6350, 0.0780, 0.1840;...
                 0, 0.4470, 0.7410; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.6350, 0.0780, 0.1840];
cfg.linewidth = 1.5;
ft_singleplotER(cfg, Ret_GA_Cat_Insect, Ret_GA_Cat_Bird, Ret_GA_Cat_Mammal, Ret_GA_Cat_Sea, Ret_GA_Cat_Fruit, Ret_GA_Cat_Electronic, Ret_GA_Cat_Veggie, Ret_GA_Cat_Clothe)
cfg.channel = ROI_temp';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Ret_GA_Cat_Insect, Ret_GA_Cat_Bird, Ret_GA_Cat_Mammal, Ret_GA_Cat_Sea, Ret_GA_Cat_Fruit, Ret_GA_Cat_Electronic, Ret_GA_Cat_Veggie, Ret_GA_Cat_Clothe)


% Difference Perceptual and Semantic

Ret_GA_Perc_Diff = Ret_GA_Perc_Drawing;
Ret_GA_Perc_Diff.avg = Ret_GA_Perc_Drawing.avg - Ret_GA_Perc_Picture.avg;

Ret_GA_Sem_Diff = Ret_GA_Sem_Animate;
Ret_GA_Sem_Diff.avg = Ret_GA_Sem_Animate.avg - Ret_GA_Sem_NonAnimate.avg;

cfg = [];
cfg.parameter = 'avg';
%cfg.baseline = [-0.2 0];
cfg.channel = ROI_occ';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Ret_GA_Perc_Diff, Ret_GA_Sem_Diff)
cfg.channel = ROI_temp';
cfg.showlegend = 'yes';
ft_singleplotER(cfg, Ret_GA_Perc_Diff, Ret_GA_Sem_Diff)




figure
plot([]); hold on
plot([Enc_Data_EEG.TimeVec1024(1) Enc_Data_EEG.TimeVec1024(end)],[0 0],'--k','linewidth',1.5)
plot(Enc_Data_EEG.TimeVec1024, squeeze(Enc_Data_EEG.Data(57,1,:)),'b','linewidth',1)
set(gca,'xlim',[Enc_Data_EEG.TimeVec1024(1) Enc_Data_EEG.TimeVec1024(end)],'ydir','reverse')
hold off; box off



