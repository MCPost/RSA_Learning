
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

%save('RSA_Data_Enc', 'Subj_names', 'elecs', 'ROI_all_idx')


c = parcluster();
jobHandles = cell(length(Subj_names),1);
for sub = 1:length(Subj_names)
    jobHandles{sub}    = batch(c, @RSA_Enc_GA_LoadData, 0, {sub,Subj_names});
end

jobHandles{22}
diary(jobHandles{22})

%for sub = 1:length(Subj_names)
%    delete(jobHandles{sub})
%end
