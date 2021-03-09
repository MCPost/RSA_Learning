%% Frequency Methods Permutations Per Sub


elecs        = Elec_Info.elecs;
ROI_all_idx  = Elec_Info.ROI_all_idx;
ROI_occ_idx  = Elec_Info.ROI_occ_idx;
ROI_temp_idx = Elec_Info.ROI_temp_idx;


freq_meth = {'MTM','Wavelet','FHilbert'};
freq_band = {'theta','alpha','beta','gamma1','gamma2'};
centerfreq = [6, 10, 20, 40, 75];
lu_filter_bound = [2 4; 8 12; 15 25; 30 50; 65 85];
tapsmofrq = [2 2 5 10 15];
BL_perms = {'NoBL','logBL','pctBL'};


%% Import Data
eeg_cfg = [];
eeg_cfg.Name = Name;
eeg_cfg.chan_label = elecs; 
eeg_cfg.chan_idx = ROI_all_idx'; 
eeg_cfg.Art_corr = 0;
eeg_cfg.BL_corr = 'no';
eeg_cfg.BL_wind = [];
Data_EEG = createEEGData(eeg_cfg);


Preproc_Data_1 = load(['Prepro_2_encoding_full_',Name,'.mat'], ['r_1_',Name]);
Preproc_Data_1 = Preproc_Data_1.(cell2mat(fieldnames(Preproc_Data_1)));
Preproc_Data_1.trial = Data_EEG.Data;

%% Zero pad to trials
zpData1 = Preproc_Data_1;
zpData1.time = Preproc_Data_1.time(1)-floor(2353/2)/1024:1/1024:Preproc_Data_1.time(end)+ceil(2353/2)/1024;
zpData1.trial = zeros(size(Preproc_Data_1.trial)+[0 0 2353]);
for tp = 1:length(Preproc_Data_1.time)
    zpData1.trial(:,:,floor(2353/2)+tp) = Preproc_Data_1.trial(:,:,tp);
end

%% Fieldtrip Analysis

cycles = 7;

rsa_cfg = Data_EEG;
rsa_cfg.slide_window = 0.080; 
rsa_cfg.slide_step = 0.004;
rsa_cfg.window_average = 'gaussian';
rsa_cfg.meas128 = '';
rsa_cfg.meas16 = {'LDA','SVM','euclidian'};
rsa_cfg.Cktl_blank_rm = false;
rsa_cfg.only16 = true;
rsa_cfg.MNN = true;
for frq = 1:length(centerfreq)
    
    cfg = [];
    cfg.method      = 'mtmconvol';
    cfg.output      = 'pow';
    cfg.pad         = 'nextpow2';
    cfg.keeptrials  = 'yes';
    cfg.foi         = centerfreq(frq);
    cfg.taper       = 'dpss';
    cfg.tapsmofrq   = tapsmofrq(frq);
    cfg.t_ftimwin   = cycles./cfg.foi;
    cfg.toi         = Preproc_Data_1.time;
    cfg.feedback    = 'no';

    [freq] = ft_freqanalysis(cfg, zpData1);
    
    %% Baseline
    
    for bl = 1:3
        if(bl == 1)      % No Baseline
            rsa_cfg.Data = squeeze(freq.powspctrm);
        elseif(bl == 2)  % log Baseline
            rsa_cfg.Data = 10*log10(bsxfun(@rdivide, squeeze(freq.powspctrm), squeeze(nanmean(freq.powspctrm(:,:,:,1:156),4))));
        else             % Percent Baseline
            rsa_cfg.Data = bsxfun(@rdivide, bsxfun(@minus, squeeze(freq.powspctrm), squeeze(nanmean(freq.powspctrm(:,:,:,1:156),4))), squeeze(nanmean(freq.powspctrm(:,:,:,1:156),4)));
        end
        % Occipital
        rsa_cfg.curROI = ROI_occ_idx;
        rsa_cfg.curROI_name = 'occipital';
        a = createRSA(rsa_cfg);
        % Temporal
        rsa_cfg.curROI = ROI_temp_idx;
        rsa_cfg.curROI_name = 'temporal';
        a = createRSA(rsa_cfg);
    end
    
end

data_hann_nobl   = squeeze(nanmean(nanmean(freq.powspctrm(:,ROI_occ_idx,:,:),2),1))'; 
data_hann_logbl  = squeeze(nanmean(nanmean(10*log10(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))),2),1))';
data_hann_prctbl = squeeze(nanmean(nanmean(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))*100,2),1))';





%% Wavelet Analysis

centerfreq =  75;

freq = [];
freq.time = Preproc_Data_1.time;
freq.freq = centerfreq;
tmp = zeros(size(Preproc_Data_1.trial,1),size(Preproc_Data_1.trial,3),size(zpData1.trial,2),length(centerfreq));
for frq = 1:length(centerfreq)

    % setup wavelet convolution and outputs
    time = -1:1/1024:1;
    half_of_wavelet_size = (length(time)-1)/2;

    % FFT parameters
    n_wavelet      = length(time);
    n_data         = size(zpData1.trial,3)*size(zpData1.trial,1);
    n_convolution  = n_wavelet+n_data-1;
    wavelet_cycles = 7;

    for chan = 1:size(zpData1.trial,2)

        % FFT of data (note: this doesn't change on frequency iteration)
        fft_data = fft(reshape(squeeze(zpData1.trial(:,chan,:))',1,n_data),n_convolution);

        % create wavelet and run convolution
        fft_wavelet             = fft(exp(2*1i*pi*centerfreq(frq).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq(frq)))^2)),n_convolution);
        convolution_result_fft  = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq(frq)));
        convolution_result_fft  = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        convolution_result_fft  = reshape(convolution_result_fft,size(zpData1.trial,3),size(zpData1.trial,1))';
        tmp(:,:,chan,frq)       = abs(convolution_result_fft(:,floor(2353/2)+(1:length(Preproc_Data_1.time)))).^2;  
    end
    freq.powspctrm  = permute(tmp,[1 3 4 2]);
end




%% Filter-Hilbert Method


nyquist = 1024/2;
lu_filter_bound = [65 85]; % Hz
transition_width   = 0.1;

tmp = zeros(length(Preproc_Data_1.time),size(zpData1.trial,1),size(zpData1.trial,2),size(lu_filter_bound,1));
for frq = 1:size(lu_filter_bound,1)
    
    filter_order       = round(3*(1024/lu_filter_bound(frq,1)));
    % create the filter shape (this is explained more in the text around figure 14.4)
    ffrequencies  = [ 0 (1-transition_width)*lu_filter_bound(frq,1) lu_filter_bound(frq,1) lu_filter_bound(frq,2) (1+transition_width)*lu_filter_bound(frq,2) nyquist ]/nyquist;
    idealresponse = [ 0 0 1 1 0 0 ];
    filterweights = firls(filter_order,ffrequencies,idealresponse);

    % apply the filter kernal to the data to obtain the band-pass filtered signal
    freq = [];
    freq.time = Preproc_Data_1.time;
    for chani=1:size(zpData1.trial,2)
        data2filter_cat = squeeze(double(reshape(squeeze(zpData1.trial(:,chani,:))',1,size(zpData1.trial,3)*size(zpData1.trial,1))));
        filtdat_cat = reshape(filtfilt(filterweights,1,data2filter_cat),size(zpData1.trial,3),size(zpData1.trial,1));
        tmp(:,:,chani,frq) = abs(hilbert(filtdat_cat(floor(2353/2)+(1:length(Preproc_Data_1.time)),:)));
    end
end
freq.powspctrm  = permute(tmp,[2 3 4 1]);

data_fhil_own_nobl   = squeeze(nanmean(nanmean(freq.powspctrm(:,ROI_occ_idx,:,:),2),1))'; 
data_fhil_own_logbl  = squeeze(nanmean(nanmean(10*log10(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))),2),1))';
data_fhil_own_prctbl = squeeze(nanmean(nanmean(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))*100,2),1))';

