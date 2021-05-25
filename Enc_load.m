
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

% ROI_idx = {ROI_occ_idx, 'Occipital';...
%            ROI_temp_idx, 'Temporal';...
%            ROI_front_idx, 'Frontal';...
%            ROI_pari_idx, 'Parietal'};
% save_name = 'RSA_Data/RSA_Data_Enc_';
% 
% c = parcluster();
% jobHandles = cell(length(Subj_names),1);
% for sub = 1:length(Subj_names)
%     save([save_name,Subj_names{sub}], 'Subj_names', 'elecs', 'ROI_idx')
%     jobHandles{sub}    = batch(c, @RSA_Enc_GA_LoadData, 0, {sub,Subj_names, ROI_idx, [save_name,Subj_names{sub}]});
% end

ROI_idx = {ROI_all_idx, 'Wholehead'};
save_name = 'RSA_Data_Enc_wholehead';

save(save_name, 'Subj_names', 'elecs', 'ROI_idx')

c = parcluster();
jobHandles = cell(length(Subj_names),1);
for sub = 1:length(Subj_names)
    jobHandles{sub}    = batch(c, @RSA_Enc_GA_LoadData, 0, {sub,Subj_names, ROI_idx, save_name});
end

jobHandles{22}
diary(jobHandles{22})

%for sub = 1:length(Subj_names)
%    delete(jobHandles{sub})
%end
