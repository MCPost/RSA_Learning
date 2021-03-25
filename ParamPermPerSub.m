%% Compute Parameter Permutations per Subject

function ParamPermPerSub(Name, Elec_Info)

elecs        = Elec_Info.elecs;
ROI_all_idx  = Elec_Info.ROI_all_idx;
ROI_occ_idx  = Elec_Info.ROI_occ_idx;
ROI_temp_idx = Elec_Info.ROI_temp_idx;


artef_perms = {'ArtCorr','NoArtCorr'};
BL_perms = {'NoBL','MinBL','CovBL'};
CBR_perms = {'CBR','NoCBR'};
MNN_perms = {'NoMNN','MNN'};


%% 1. Permutation: Artifact Rejection (yes / no)

for art = 1:length(artef_perms)
    
    art_corr = logical(mod(art,2));
    
%% 2. Permutation: Baseline Correction (no / minus BL / covariate BL)

    for blc = 1:length(BL_perms)
        
        bl_corr = 'no';
        if(blc == 2)
            bl_corr = 'minBL';
        elseif(blc == 3)
            bl_corr = 'covBL';
        end


        %% Import Data
        eeg_cfg = [];
        eeg_cfg.Name = Name;
        eeg_cfg.chan_label = elecs; 
        eeg_cfg.chan_idx = ROI_all_idx'; 
        eeg_cfg.Art_corr = art_corr;
        eeg_cfg.BL_corr = bl_corr;
        eeg_cfg.BL_wind = [-0.26 -0.05];
        Data_EEG.(['Data_EEG_',artef_perms{art},'_',BL_perms{blc}]) = createEEGData(eeg_cfg);

        save(['RSA_DiffParams/Data_EEG_',Name,'_',artef_perms{art},'_',BL_perms{blc}],'-struct','Data_EEG','-append')
        Data_EEG = Data_EEG.(['Data_EEG_',artef_perms{art},'_',BL_perms{blc}]);


%% 3. Permutation: Cocktail Blank Removal (yes / no)

        for cbr = 1:length(CBR_perms)

            cocktail_blank = logical(mod(cbr,2));

            fprintf('\n')
            fprintf('  %s (%i/%i) --- %s (%i/%i) --- %s (%i/%i) --- Subject: %s started!',artef_perms{art}, art, length(artef_perms), BL_perms{blc}, blc, length(BL_perms), CBR_perms{cbr}, cbr, length(CBR_perms), Name)
            fprintf('\n')
            

            %% RSA Function
            rsa_cfg = Data_EEG;
            rsa_cfg.slide_window = 0.080; 
            rsa_cfg.slide_step = 0.004;
            rsa_cfg.window_average = 'gaussian';
            rsa_cfg.Cktl_blank_rm = cocktail_blank;
            rsa_cfg.noMDS = false;

            rsa_cfg.only16 = false;
            rsa_cfg.MNN = false;
            rsa_cfg.curROI = ROI_occ_idx;
            rsa_cfg.curROI_name = 'occipital';
            RSA_Data_OCC.(['RSA_Data_OCC_',artef_perms{art},'_',BL_perms{blc},'_',CBR_perms{cbr},'_',MNN_perms{1}]) = createRSA(rsa_cfg);

            save(['RSA_DiffParams/RSA_Data_',Name,'_',artef_perms{art},'_',BL_perms{blc},'_',CBR_perms{cbr},'_',MNN_perms{1}],'-struct','RSA_Data_OCC','-append')
            clear('RSA_Data_OCC')

            rsa_cfg.only16 = true;
            rsa_cfg.MNN = true;
            RSA_Data_OCC.(['RSA_Data_OCC_',artef_perms{art},'_',BL_perms{blc},'_',CBR_perms{cbr},'_',MNN_perms{2}]) = createRSA(rsa_cfg);

            save(['RSA_DiffParams/RSA_Data_',Name,'_',artef_perms{art},'_',BL_perms{blc},'_',CBR_perms{cbr},'_',MNN_perms{2}],'-struct','RSA_Data_OCC','-append')
            clear('RSA_Data_OCC')

            rsa_cfg.only16 = false;
            rsa_cfg.MNN = false;
            rsa_cfg.curROI = ROI_temp_idx;
            rsa_cfg.curROI_name = 'temporal';
            RSA_Data_TMP.(['RSA_Data_TMP_',artef_perms{art},'_',BL_perms{blc},'_',CBR_perms{cbr},'_',MNN_perms{1}]) = createRSA(rsa_cfg);

            save(['RSA_DiffParams/RSA_Data_',Name,'_',artef_perms{art},'_',BL_perms{blc},'_',CBR_perms{cbr},'_',MNN_perms{1}],'-struct','RSA_Data_TMP','-append')
            clear('RSA_Data_TMP')


            rsa_cfg.only16 = true;
            rsa_cfg.MNN = true;
            RSA_Data_TMP.(['RSA_Data_TMP_',artef_perms{art},'_',BL_perms{blc},'_',CBR_perms{cbr},'_',MNN_perms{2}]) = createRSA(rsa_cfg);
            
            save(['RSA_DiffParams/RSA_Data_',Name,'_',artef_perms{art},'_',BL_perms{blc},'_',CBR_perms{cbr},'_',MNN_perms{2}],'-struct','RSA_Data_TMP','-append')
            clear('RSA_Data_TMP')

        end
        clear Data_EEG    
                

    end
    

end

end