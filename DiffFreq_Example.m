%% Different Time-Frequency splits Example

%% Subject Names
Subj_names = {'AG','CEWD','CM','DS','FVM'};


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
elecs = {'A1';'A2';'A3';'A4';'A5';'A6';'A7';'A8';'A9';'A10';'A11';'A12';'A13';'A14';'A15';'A16';'A17';...
          'A18';'A19';'A20';'A21';'A22';'A23';'A24';'A25';'A26';'A27';'A28';'A29';'A30';'A31';'A32';'B1';...
          'B2';'B3';'B4';'B5';'B6';'B7';'B8';'B9';'B10';'B11';'B12';'B13';'B14';'B15';'B16';'B17';'B18';...
          'B19';'B20';'B21';'B22';'B23';'B24';'B25';'B26';'B27';'B28';'B29';'B30';'B31';'B32';'C1';'C2';...
          'C3';'C4';'C5';'C6';'C7';'C8';'C9';'C10';'C11';'C12';'C13';'C14';'C15';'C16';'C17';'C18';'C19';...
          'C20';'C21';'C22';'C23';'C24';'C25';'C26';'C27';'C28';'C29';'C30';'C31';'C32';'D1';'D2';'D3';...
          'D4';'D5';'D6';'D7';'D8';'D9';'D10';'D11';'D12';'D13';'D14';'D15';'D16';'D17';'D18';'D19';'D20';...
          'D21';'D22';'D23';'D24';'D25';'D26';'D27';'D28';'D29';'D30';'D31';'D32'};

% All Electrodes
ROI_all_idx = 1:128;

% Occipital Cortex ROI
ROI_occ = {'A8','A9','A10','A11','A12','A13','A14','A15','A16','A17','A21','A22','A23','A24','A25','A26',...
           'A27','A28','A29','A30','B5','B6','B7','B8','B9',...
           'D32','D31','D30','D29','A6','A7','A5','A18','A19','A20','A32','A31','B3','B4','B13','B12','B11','B10'};
ROI_occ_idx = find(cell2mat(cellfun(@(x) any(strcmp(x, ROI_occ)), elecs, 'UniformOutput', 0)));

% Temporal Cortex ROI
ROI_temp = {'B18','B17','B16','B15','B14','B22','B23','B24','B25','B26','B31','B30','B29','B28','B27','C3',...
            'C4','C5','C6','C7','D28','D27','D26','D25','D24','D19','D20','D21','D22','D23','D8','D9','D10',...
            'D11','D12','D7','D6','D5','D4','D3'};
ROI_temp_idx = find(cell2mat(cellfun(@(x) any(strcmp(x, ROI_temp)), elecs, 'UniformOutput', 0)));

Name = 'AG';

Preproc_Data_1 = load(['Prepro_2_encoding_full_',Name,'.mat'], ['r_1_',Name]);
Preproc_Data_1 = Preproc_Data_1.(cell2mat(fieldnames(Preproc_Data_1)));
Preproc_Data_2 = load(['Prepro_2_encoding_full_',Name,'.mat'], ['r_2_',Name]);
Preproc_Data_2 = Preproc_Data_2.(cell2mat(fieldnames(Preproc_Data_2)));

timeVec = Preproc_Data_1.time;


%% Zero pad to trials
zpData1 = Preproc_Data_1;
zpData1.time = Preproc_Data_1.time(1)-floor(2353/2)/1024:1/1024:Preproc_Data_1.time(end)+ceil(2353/2)/1024;
zpData1.trial = zeros(size(Preproc_Data_1.trial)+[0 0 2353]);
zpData2 = Preproc_Data_2;
zpData2.time = Preproc_Data_2.time(1)-floor(2353/2)/1024:Preproc_Data_2.time(end)+ceil(2353/2)/1024;
zpData2.trial = zeros(size(Preproc_Data_2.trial)+[0 0 2353]);
for tp = 1:length(Preproc_Data_1.time)
    zpData1.trial(:,:,floor(2353/2)+tp) = Preproc_Data_1.trial(:,:,tp);
    zpData2.trial(:,:,floor(2353/2)+tp) = Preproc_Data_2.trial(:,:,tp);
end
%squeeze(Preproc_Data_1.trial(1,1,:))
%squeeze(zpdat(1,1,:))


%% Fieldtrip Analysis

cycles = 7;

cfg = [];
cfg.method      = 'mtmconvol';
cfg.output      = 'pow';
cfg.pad         = 'nextpow2';
cfg.keeptrials  = 'yes';
cfg.foi         = 75;%[2 3 4 5 6 7 8 9 10 12 14 16 18 20 25 30 35 40 50 60 70 80 90];
cfg.taper       = 'dpss' ;
cfg.tapsmofrq   = 20;
cfg.t_ftimwin   = cycles./cfg.foi;
cfg.toi         = Preproc_Data_1.time; 

[freq] = ft_freqanalysis(cfg, zpData1);


data_hann_nobl   = squeeze(nanmean(nanmean(freq.powspctrm(:,ROI_occ_idx,:,:),2),1))'; 
data_hann_logbl  = squeeze(nanmean(nanmean(10*log10(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))),2),1))';
data_hann_prctbl = squeeze(nanmean(nanmean(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))*100,2),1))';


figure
subplot(3,1,1)
plot(freq.time, data_hann_nobl)
set(gca,'xlim',[timeVec(1) timeVec(end)])
subplot(3,1,2)
plot(freq.time, data_hann_logbl)
set(gca,'xlim',[timeVec(1) timeVec(end)])
subplot(3,1,3)
plot(freq.time, data_hann_prctbl)
set(gca,'xlim',[timeVec(1) timeVec(end)])


% figure
% subplot(3,1,1:2)
% plot_data = squeeze(nanmean(nanmean(freq.powspctrm(:,ROI_occ_idx,:,:),1),2));
% contourf(freq.time,freq.freq,10*log10(bsxfun(@rdivide,plot_data, nanmean(plot_data(:,1:206),2))),40,'linestyle','none')
% set(gca,'yscale','log','ytick',cfg.foi,'Fontsize',9)
% subplot(3,1,3)
% plot(freq.time, squeeze(nanmean(plot_data,1)))
% set(gca,'xlim',[freq.time(1) freq.time(end)])
% 
% figure
% plot_data = squeeze(nanmean(nanmean(freq.powspctrm(:,ROI_occ_idx,:,:),1),2))';
% plot(freq.time, plot_data)
% set(gca,'xlim',[freq.time(1) freq.time(end)])



%% Wavelet Analysis

centerfreq =  75; %[2 3 4 5 6 7 8 9 10 12 14 16 18 20 25 30 35 40 50 60 70 80 90]; % in Hz

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
end
freq.powspctrm  = permute(tmp,[1 3 4 2]);

data_wvl_own_nobl   = squeeze(nanmean(nanmean(freq.powspctrm(:,ROI_occ_idx,:,:),2),1))'; 
data_wvl_own_logbl  = squeeze(nanmean(nanmean(10*log10(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))),2),1))';
data_wvl_own_prctbl = squeeze(nanmean(nanmean(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))*100,2),1))';

figure
subplot(3,2,1)
plot(freq.time, data_wvl_own_nobl)
set(gca,'xlim',[timeVec(1) timeVec(end)])
subplot(3,2,3)
plot(freq.time, data_wvl_own_logbl)
set(gca,'xlim',[timeVec(1) timeVec(end)])
subplot(3,2,5)
plot(freq.time, data_wvl_own_prctbl)
set(gca,'xlim',[timeVec(1) timeVec(end)])


cfg = [];
cfg.method      = 'wavelet';
cfg.output      = 'pow';
cfg.pad         = 'nextpow2';
cfg.keeptrials  = 'yes';
cfg.foi         = 75;%[2 3 4 5 6 7 8 9 10 12 14 16 18 20 25 30 35 40 50 60 70 80 90];
cfg.width       = cycles;
cfg.toi         = Preproc_Data_1.time; 

[freq] = ft_freqanalysis(cfg, zpData1);

data_wvl_ft_nobl   = squeeze(nanmean(nanmean(freq.powspctrm(:,ROI_occ_idx,:,:),2),1))'; 
data_wvl_ft_logbl  = squeeze(nanmean(nanmean(10*log10(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))),2),1))';
data_wvl_ft_prctbl = squeeze(nanmean(nanmean(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))*100,2),1))';

subplot(3,2,2)
plot(freq.time, data_wvl_ft_nobl)
set(gca,'xlim',[timeVec(1) timeVec(end)])
subplot(3,2,4)
plot(freq.time, data_wvl_ft_logbl)
set(gca,'xlim',[timeVec(1) timeVec(end)])
subplot(3,2,6)
plot(freq.time, data_wvl_ft_prctbl)
set(gca,'xlim',[timeVec(1) timeVec(end)])

% figure
% subplot(3,1,1:2)
% plot_data = squeeze(nanmean(nanmean(freq.powspctrm(:,ROI_occ_idx,:,:),1),2));
% contourf(freq.time,freq.freq,10*log10(bsxfun(@rdivide,plot_data, nanmean(plot_data(:,1:206),2))),40,'linestyle','none'); colorbar
% set(gca,'yscale','log','ytick',[2 4 6 8 10 12 15 20 25 30 50 70],'Fontsize',9)
% subplot(3,1,3)
% plot(freq.time, squeeze(nanmean(plot_data,1)))
% set(gca,'xlim',[freq.time(1) freq.time(end)])
% 
% figure
% plot_data = squeeze(nanmean(nanmean(freq.powspctrm(:,ROI_occ_idx,:,:),1),2))';
% plot(freq.time, plot_data)
% set(gca,'xlim',[freq.time(1) freq.time(end)])

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

figure
subplot(3,2,1)
plot(freq.time, data_fhil_own_nobl)
set(gca,'xlim',[timeVec(1) timeVec(end)])
subplot(3,2,3)
plot(freq.time, data_fhil_own_logbl)
set(gca,'xlim',[timeVec(1) timeVec(end)])
subplot(3,2,5)
plot(freq.time, data_fhil_own_prctbl)
set(gca,'xlim',[timeVec(1) timeVec(end)])

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfilttype = 'fir';
tmp = zeros(size(Preproc_Data_1.trial,2),size(Preproc_Data_1.trial,3),size(Preproc_Data_1.trial,1),size(lu_filter_bound,1));
for frq = 1:size(lu_filter_bound,1)
    cfg.bpfreq = lu_filter_bound(frq,:);
    [data_bp] = ft_preprocessing(cfg, Preproc_Data_1);
    for trl = 1:size(data_bp.trial,1)
        tmp(:,:,trl,frq) = ft_preproc_hilbert(squeeze(data_bp.trial(trl,:,:)));
    end
end
freq.powspctrm  = permute(tmp,[3 1 4 2]);

data_fhil_ft_nobl   = squeeze(nanmean(nanmean(dat(:,ROI_occ_idx,:),2),1))'; 
data_fhil_ft_logbl  = squeeze(nanmean(nanmean(10*log10(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))),2),1))';
data_fhil_ft_prctbl = squeeze(nanmean(nanmean(bsxfun(@rdivide,squeeze(freq.powspctrm(:,ROI_occ_idx,:,:)), squeeze(nanmean(freq.powspctrm(:,ROI_occ_idx,:,1:156),4)))*100,2),1))';


subplot(3,2,2)
plot(freq.time, data_fhil_ft_nobl)
set(gca,'xlim',[timeVec(1) timeVec(end)],'ylim',[1 1.8])
subplot(3,2,4)
plot(freq.time, data_fhil_ft_logbl)
set(gca,'xlim',[timeVec(1) timeVec(end)],'ylim',[-2 0.5])
subplot(3,2,6)
plot(freq.time, data_fhil_ft_prctbl)
set(gca,'xlim',[timeVec(1) timeVec(end)],'ylim',[80 130])

% figure
% plot_data = squeeze(nanmean(nanmean(freq.powspctrm(:,ROI_occ_idx,:),1),2))';
% plot(freq.time, plot_data)
% set(gca,'xlim',[freq.time(1) freq.time(end)])


