%% Single Subject try out retrieval

%% Subject Names
[~,message,~] = fileattrib('Preproc_EEG_Data\Retrieval_response_locked\*');
currentdir = pwd;
filenames = strrep({message([message.directory] == 0).Name}',[currentdir,'\Preproc_EEG_Data\Retrieval_response_locked\'],'');
Subj_names = cellfun(@(x) x{1}(10:end), regexp(filenames,'_(\w*).','tokens','once'),'UniformOutput', 0);

sub = 10;

%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs


%% Import Data
load(['Preproc_EEG_Data/Retrieval_response_locked/Ret_EEG_Data_',Subj_names{sub}])


% Plot Grand Averages
figure
plot(Ret_Data_EEG.TimeVec1024, squeeze(nanmean(nanmean(Ret_Data_EEG.Data(:,ROI_occ_idx,:),2),1)),'b','linewidth',2)
hold on
plot(Ret_Data_EEG.TimeVec1024, squeeze(nanmean(nanmean(Ret_Data_EEG.Data(:,ROI_temp_idx,:),2),1)),'r','linewidth',2)
set(gca,'linewidth',1.5,'xlim',[-4 0.2]); box off
hold off

% Conditions
figure
subplot(2,2,1)
plot(Ret_Data_EEG.TimeVec1024, squeeze(nanmean(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,5)) == 1,ROI_occ_idx,:),2),1)),'b','linewidth',2)
hold on
plot(Ret_Data_EEG.TimeVec1024, squeeze(nanmean(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,5)) == 2,ROI_occ_idx,:),2),1)),'r','linewidth',2)
set(gca,'linewidth',1.5,'xlim',[-4 0.2]); box off
hold off
subplot(2,2,2)
plot(Ret_Data_EEG.TimeVec1024, squeeze(nanmean(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,7)) == 1,ROI_occ_idx,:),2),1)),'b','linewidth',2)
hold on
plot(Ret_Data_EEG.TimeVec1024, squeeze(nanmean(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,7)) == 2,ROI_occ_idx,:),2),1)),'r','linewidth',2)
set(gca,'linewidth',1.5,'xlim',[-4 0.2]); box off
hold off
subplot(2,2,3)
plot(Ret_Data_EEG.TimeVec1024, squeeze(nanmean(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,5)) == 1,ROI_temp_idx,:),2),1)),'b','linewidth',2)
hold on
plot(Ret_Data_EEG.TimeVec1024, squeeze(nanmean(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,5)) == 2,ROI_temp_idx,:),2),1)),'r','linewidth',2)
set(gca,'linewidth',1.5,'xlim',[-4 0.2]); box off
hold off
subplot(2,2,4)
plot(Ret_Data_EEG.TimeVec1024, squeeze(nanmean(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,7)) == 1,ROI_temp_idx,:),2),1)),'b','linewidth',2)
hold on
plot(Ret_Data_EEG.TimeVec1024, squeeze(nanmean(nanmean(Ret_Data_EEG.Data(cell2mat(Ret_Data_EEG.TrialInfo(:,7)) == 2,ROI_temp_idx,:),2),1)),'r','linewidth',2)
set(gca,'linewidth',1.5,'xlim',[-4 0.2]); box off
hold off

