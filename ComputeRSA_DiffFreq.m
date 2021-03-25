%% Try RSA with different frequency methods


%% Subject Names
Subj_names = {'AG','CEWD','CM','DS','FVM'};


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs

% c = parcluster();
% j_AG    = batch(c, @FreqMethPerSub, 0, {Subj_names{1}, Elec_Info});
% j_CEWD  = batch(c, @FreqMethPerSub, 0, {Subj_names{2}, Elec_Info});
% j_CM    = batch(c, @FreqMethPerSub, 0, {Subj_names{3}, Elec_Info});
% j_DS    = batch(c, @FreqMethPerSub, 0, {Subj_names{4}, Elec_Info});
% j_FVM   = batch(c, @FreqMethPerSub, 0, {Subj_names{5}, Elec_Info});

%FreqMethPerSub(Subj_names{1}, Elec_Info)

diary(j_AG)
diary(j_CEWD)
diary(j_CM)
diary(j_DS)
diary(j_FVM)

%delete(j_AG)
%delete(j_CEWD)
%delete(j_CM)
%delete(j_DS)
%delete(j_FVM)



