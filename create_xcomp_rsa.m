%% Compute Crosscomparison of RSA timecourse

function CrossComp_RSA = create_xcomp_rsa(cfg, Data1, Data2)

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
Hyp_perceptual = cfg.Hyp_perceptual;
Hyp_semantic = cfg.Hyp_semantic;
TimeVec_Data1 = Data1.TimeVec;
TimeVec_Data2 = Data2.TimeVec; 
curROI = cfg.ROI;
only16 = cfg.only16;


%% Create RSA Matrix
samplingrate1 = unique(1./diff(TimeVec_Data1)); %Hz
slide_window_idx1 = round(samplingrate1*(slide_window_s/2));
slide_step_idx1 = round(samplingrate1*slide_step_s);
TimeVec1_idx = 1:slide_step_idx1:length(TimeVec_Data1);
TimeVec1 = TimeVec_Data1(TimeVec1_idx);
cfg.TimeVec1 = TimeVec1;

samplingrate2 = unique(1./diff(TimeVec_Data2)); %Hz
slide_window_idx2 = round(samplingrate2*(slide_window_s/2));
slide_step_idx2 = round(samplingrate2*slide_step_s);
TimeVec2_idx = 1:slide_step_idx2:length(TimeVec_Data2);
TimeVec2 = TimeVec_Data2(TimeVec2_idx);
cfg.TimeVec2 = TimeVec2;


% Save Data
CrossComp_RSA = cfg;
CrossComp_RSA.RSA_full   = []; 

fprintf('\n')
fprintf('Create Crosscompare RSA Matrix (Slide Window: %1.3f s | Slide Step: %1.3f s)!', slide_window_s, slide_step_s)
fprintf('\n')

if(~only16)
    
    per_ind_wi = Hyp_perceptual{1}(:) > 0;
    per_ind_bt = Hyp_perceptual{1}(:) < 0;
    sem_ind_wi = Hyp_semantic{1}(:) > 0;
    sem_ind_bt = Hyp_semantic{1}(:) < 0;
    
    for r = 1:length(curROI)
        
        CrossComp_RSA.RSA_full.(curROI{r}) = zeros(size(Data1.(curROI{r}).full_Data,1),4,length(TimeVec2),length(TimeVec1));
        fprintf('\n')
        nbytes = fprintf('ROI: %s  --  Progress Cross Correlation:  0.0 %%',curROI{r});
        for tp1 = 1:length(TimeVec1)
            
            time_window1 = max(TimeVec1_idx(tp1) - slide_window_idx1, 1):min(TimeVec1_idx(tp1) + slide_window_idx1, length(TimeVec_Data1));
            
            for tp2 = 1:length(TimeVec2)
                
                time_window2 = max(TimeVec2_idx(tp2) - slide_window_idx2, 1):min(TimeVec2_idx(tp2) + slide_window_idx2, length(TimeVec_Data2));

                CrossComp_RSA.RSA_full.(curROI{r})(:,1,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).full_Data(:,time_window1,per_ind_wi),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).full_Data(:,time_window2,per_ind_wi),2,length(time_window2)))', 'type', 'spearman'));
                CrossComp_RSA.RSA_full.(curROI{r})(:,2,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).full_Data(:,time_window1,per_ind_bt),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).full_Data(:,time_window2,per_ind_bt),2,length(time_window2)))', 'type', 'spearman'));
                CrossComp_RSA.RSA_full.(curROI{r})(:,3,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).full_Data(:,time_window1,sem_ind_wi),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).full_Data(:,time_window2,sem_ind_wi),2,length(time_window2)))', 'type', 'spearman'));
                CrossComp_RSA.RSA_full.(curROI{r})(:,4,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).full_Data(:,time_window1,sem_ind_bt),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).full_Data(:,time_window2,sem_ind_bt),2,length(time_window2)))', 'type', 'spearman'));
            
            end
            fprintf(repmat('\b',1,nbytes))
            nbytes = fprintf('ROI: %s  --  Progress Cross Correlation: %3.2f %%',curROI{r},(tp1 / length(TimeVec1))*100);
        end
        fprintf('\n') 
    end
    
end


per_ind_wi = Hyp_perceptual{2}(:) > 0;
per_ind_bt = Hyp_perceptual{2}(:) < 0;
sem_ind_wi = Hyp_semantic{2}(:) > 0;
sem_ind_bt = Hyp_semantic{2}(:) < 0;

CrossComp_RSA.RSA_red16   = []; 

for r = 1:length(curROI)

    CrossComp_RSA.RSA_red16.(curROI{r}) = zeros(size(Data1.(curROI{r}).red16_Data,1),4,length(TimeVec2),length(TimeVec1));
    fprintf('\n')
    nbytes = fprintf('ROI: %s  --  Progress Cross Correlation:  0.0 %%',curROI{r});
    for tp1 = 1:length(TimeVec1)
        
        time_window1 = max(TimeVec1_idx(tp1) - slide_window_idx1, 1):min(TimeVec1_idx(tp1) + slide_window_idx1, length(TimeVec_Data1));
        
        for tp2 = 1:length(TimeVec2)
            
            time_window2 = max(TimeVec2_idx(tp2) - slide_window_idx2, 1):min(TimeVec2_idx(tp2) + slide_window_idx2, length(TimeVec_Data2));
            
            CrossComp_RSA.RSA_red16.(curROI{r})(:,1,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).red16_Data(:,time_window1,per_ind_wi),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,per_ind_wi),2,length(time_window2)))', 'type', 'spearman'));
            CrossComp_RSA.RSA_red16.(curROI{r})(:,2,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).red16_Data(:,time_window1,per_ind_bt),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,per_ind_bt),2,length(time_window2)))', 'type', 'spearman'));
            CrossComp_RSA.RSA_red16.(curROI{r})(:,3,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).red16_Data(:,time_window1,sem_ind_wi),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,sem_ind_wi),2,length(time_window2)))', 'type', 'spearman'));
            CrossComp_RSA.RSA_red16.(curROI{r})(:,4,tp2,tp1) = diag(corr(squeeze(average_kern(Data1.(curROI{r}).red16_Data(:,time_window1,sem_ind_bt),2,length(time_window1)))', squeeze(average_kern(Data2.(curROI{r}).red16_Data(:,time_window2,sem_ind_bt),2,length(time_window2)))', 'type', 'spearman'));

            %CrossComp_RSA_Mat(tp_ret,tp_enc) = fast_corr(squeeze(Data1.(ROI{1}).red16_Data(sub,tp_enc,mat_ind)), squeeze(Data2.(ROI{1}).red16_Data(sub,tp_ret,mat_ind)));
            %CrossComp_RSA_Mat(1,tp_ret,tp_enc) = corr(squeeze(Data1.(ROI{2}).red16_Data(sub,tp_enc,mat_ind_wi)), squeeze(Data2.(ROI{2}).red16_Data(sub,tp_ret,mat_ind_wi)), 'type', 'spearman');
            %CrossComp_RSA_Mat(2,tp_ret,tp_enc) = corr(squeeze(Data1.(ROI{2}).red16_Data(sub,tp_enc,mat_ind_bt)), squeeze(Data2.(ROI{2}).red16_Data(sub,tp_ret,mat_ind_bt)), 'type', 'spearman');
            
        end
        fprintf(repmat('\b',1,nbytes))
        nbytes = fprintf('ROI: %s  --  Progress Cross Correlation: %3.2f %%',curROI{r},(tp1 / length(TimeVec1))*100);
    end
    fprintf('\n')
end


%% Save Data
CrossComp_RSA.xcomp_dim        = 'sub_hyp_rettime_enctime';


end
