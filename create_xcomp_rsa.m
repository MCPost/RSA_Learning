%% Compute Crosscomparison of RSA timecourse

function CrossComp_RSA = create_xcomp_rsa(cfg)

% Parameter Settings
slide_window_s = cfg.slide_window; %s
slide_step_s = cfg.slide_step; %s
if(strcmp(cfg.window_average,'gaussian'))
    average_kern = @(x, dim, gsize) sum(bsxfun(@times, x, (gausswin(gsize)./sum(gausswin(gsize)))'),dim);
    av = 'gaussian'; % Gaussian Average
else
    average_kern = @(x, dim, gsize) mean(x,dim);
    av = 'uniform';  % Uniform Average
end
Name = cfg.Name;
TrialInfo = cfg.TrialInfo;
TimeVec1024 = cfg.TimeVec1024; 
curROI = cfg.curROI;
meas128 = cfg.meas128;
meas16  = cfg.meas16;
MNN = cfg.MNN; 
Cktl_blank_rm = cfg.Cktl_blank_rm;
only16 = cfg.only16;
Data_EEG = cfg.Data;
cfg = rmfield(cfg,'Data');
cfg = rmfield(cfg,'dim');


%% Create RSA Matrix
samplingrate = 1024; %Hz
slide_window_idx = round(samplingrate*(slide_window_s/2));
slide_step_idx = round(samplingrate*slide_step_s);
TimeVec_idx = 1:slide_step_idx:length(TimeVec1024);
TimeVec = TimeVec1024(TimeVec_idx);
cfg.TimeVec = TimeVec;


% Save Data
CrossComp_RSA = cfg;
CrossComp_RSA.RSA_full   = []; 


if(~only16)

    fprintf('\n')
    fprintf('Create Crosscompare RSA Matrix (Slide Window: %1.3f s | Slide Step: %1.3f s) with %s window averaging for Subject %s!', slide_window_s, slide_step_s, av, Name)
    fprintf('\n')
    
    
end


end
