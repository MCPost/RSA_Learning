%% Try RSA with different parameters

%% Subject Names
Subj_names = {'AG','CEWD','CM','DS','FVM'};


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs


artef_perms = {'ArtCorr','NoArtCorr'};
BL_perms = {'NoBL','MinBL','CovBL'};
CBR_perms = {'CBR','NoCBR'};
MNN_perms = {'MNN','NoMNN'};

% % Initialize Mat files
% save('RSA_DiffParams/Elec_Info.mat','Elec_Info')
% for sub = 1:length(Subj_names)
%     Name = Subj_names{sub};
%     
%     for art = 2:length(artef_perms)
%     
%         art_corr = logical(mod(art,2));
% 
%         for blc = 1:length(BL_perms)
% 
%             bl_corr = 'no';
%             if(blc == 2)
%                 bl_corr = 'minBL';
%             elseif(blc == 3)
%                 bl_corr = 'covBL';
%             end
%     
%             save(['RSA_DiffParams/Data_EEG_',Name,'_',artef_perms{art},'_',BL_perms{blc}], 'Name')
%             
%             for cbr = 1:length(CBR_perms)
% 
%                 save(['RSA_DiffParams/RSA_Data_',Name,'_',artef_perms{art},'_',BL_perms{blc},'_',CBR_perms{cbr},'_',MNN_perms{1}], 'Name')
%                 save(['RSA_DiffParams/RSA_Data_',Name,'_',artef_perms{art},'_',BL_perms{blc},'_',CBR_perms{cbr},'_',MNN_perms{2}], 'Name')
%                 
%             end
%         end
%     end
% end


%parfor sub = 1:length(Subj_names)
%    
%    ParamPermPerSub(Subj_names{sub}, Elec_Info)
%    
%end

%c = parcluster();
%j_AG    = batch(c, @ParamPermPerSub, 0, {Subj_names{1}, Elec_Info});
%j_CEWD  = batch(c, @ParamPermPerSub, 0, {Subj_names{2}, Elec_Info});
%j_CM    = batch(c, @ParamPermPerSub, 0, {Subj_names{3}, Elec_Info});
%j_DS    = batch(c, @ParamPermPerSub, 0, {Subj_names{4}, Elec_Info});
%j_FVM   = batch(c, @ParamPermPerSub, 0, {Subj_names{5}, Elec_Info});

%diary(j_AG)
%diary(j_CEWD)
%diary(j_CM)
%diary(j_DS)
%diary(j_FVM)
