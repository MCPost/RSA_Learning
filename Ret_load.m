
%% Subject Names
[~,message,~] = fileattrib('Preproc_EEG_Data/Retrieval_response_locked/*');
currentdir = pwd;
filenames = strrep({message([message.directory] == 0).Name}',[currentdir,'\Preproc_EEG_Data\Retrieval_response_locked\'],'');
Subj_names = cellfun(@(x) x{1}(10:end), regexp(filenames,'_(\w*).mat','tokens','once'),'UniformOutput', 0);


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs


save('RSA_Data_Ret_wholehead', 'Subj_names', 'elecs', 'ROI_all_idx')

c = parcluster();
jobHandles = cell(length(Subj_names),1);
for sub = 1:length(Subj_names)
    jobHandles{sub}    = batch(c, @RSA_Ret_GA_LoadData, 0, {sub, Subj_names});
end

jobHandles{22}
diary(jobHandles{22})

%for sub = 1:length(Subj_names)
%    delete(jobHandles{sub})
%end


