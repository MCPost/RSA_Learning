%% RSA with time-frequency data

%% Subject Names
Subj_names = {'AG','CEWD','CM','DS','FVM'};


%% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs


%% Load Subject Data

% ROI = {'OCC','TMP'};
% frqmth = {'FT','WL','FH'};
% frqs = {'theta','alpha','beta','gamma1','gamma2'};
% blnames = {'NoBL','logBL','pctBL'};
% measures = {'LDA','SVM','euclidian'};
% 
% for sub = 1:length(Subj_names)
%     for fm = frqmth
%         curData = load(['RSA_DiffFreq/RSA_Power_',Subj_names{sub}], ['RSA_Power_',fm{1}]);
%         for fq = frqs
%             for bl = blnames
%                 for m = 1:length(measures)
%                     if(sub == 1)
%                         dims = size(curData.(['RSA_Power_',fm{1}]).(ROI{1}).(fq{1}).(bl{1}).RSA_16{1,m});
%                         RSA_DiffFreq.(ROI{1}).(fm{1}).(fq{1}).(bl{1}){1,m} = nan(length(Subj_names), dims(3), dims(1), dims(2));
%                         RSA_DiffFreq.(ROI{2}).(fm{1}).(fq{1}).(bl{1}){1,m} = nan(length(Subj_names), dims(3), dims(1), dims(2));
%                     end
%                     RSA_DiffFreq.(ROI{1}).(fm{1}).(fq{1}).(bl{1}){1,m}(sub,:,:,:) = permute(curData.(['RSA_Power_',fm{1}]).(ROI{1}).(fq{1}).(bl{1}).RSA_16{1,m},[3 1 2]);
%                     RSA_DiffFreq.(ROI{2}).(fm{1}).(fq{1}).(bl{1}){1,m}(sub,:,:,:) = permute(curData.(['RSA_Power_',fm{1}]).(ROI{2}).(fq{1}).(bl{1}).RSA_16{1,m},[3 1 2]);
%                 end
%             end
%         end
%     end
% end
% RSA_DiffFreq.TimeVec = curData.(['RSA_Power_',fm{1}]).(ROI{1}).(fq{1}).(bl{1}).TimeVec;
% 
% 
% % Create Hypothesis Matrix
% trl_mat = [kron([1;2],ones(64,1)) kron([1;2;1;2],ones(32,1))];
% RSA_DiffFreq.Prcp_Hyp128 = zeros(size(trl_mat,1));
% RSA_DiffFreq.Sem_Hyp128 = zeros(size(trl_mat,1));
% for i = 1:size(trl_mat,1)-1
%     for j = (i+1):size(trl_mat,1)
%         if(j ~= size(trl_mat,1) - (i - 1))
%             if(trl_mat(i,1) == 1 && trl_mat(j,1) == 1)
%                 RSA_DiffFreq.Prcp_Hyp128(i,j) = 1;
%             elseif(trl_mat(i,1) == 2 && trl_mat(j,1) == 2)
%                 RSA_DiffFreq.Prcp_Hyp128(i,j) = 2;
%             else
%                 RSA_DiffFreq.Prcp_Hyp128(i,j) = -1;
%             end
% 
%             if(trl_mat(i,2) == 1 && trl_mat(j,2) == 1)
%                 RSA_DiffFreq.Sem_Hyp128(i,j) = 1;
%             elseif(trl_mat(i,2) == 2 && trl_mat(j,2) == 2)
%                 RSA_DiffFreq.Sem_Hyp128(i,j) = 2;
%             else
%                 RSA_DiffFreq.Sem_Hyp128(i,j) = -1;
%             end
%         end
%     end
% end
% 
% RSA_DiffFreq.Prcp_Hyp16 = zeros(16);
% RSA_DiffFreq.Sem_Hyp16 = zeros(16);
% for i = 1:16-1
%     for j = (i+1):16
%         if(j ~= 16 - (i - 1))
%             if(trl_mat(8*i,1) == 1 && trl_mat(8*j,1) == 1)
%                 RSA_DiffFreq.Prcp_Hyp16(i,j) = 1;
%             elseif(trl_mat(8*i,1) == 2 && trl_mat(8*j,1) == 2)
%                 RSA_DiffFreq.Prcp_Hyp16(i,j) = 2;
%             else
%                 RSA_DiffFreq.Prcp_Hyp16(i,j) = -1;
%             end
% 
%             if(trl_mat(8*i,2) == 1 && trl_mat(8*j,2) == 1)
%                 RSA_DiffFreq.Sem_Hyp16(i,j) = 1;
%             elseif(trl_mat(8*i,2) == 2 && trl_mat(8*j,2) == 2)
%                 RSA_DiffFreq.Sem_Hyp16(i,j) = 2;
%             else
%                 RSA_DiffFreq.Sem_Hyp16(i,j) = -1;
%             end
%         end
%     end
% end
% clear trl_mat
% 
% save('RSA_DiffFreq','RSA_DiffFreq','Subj_names','ROI_occ','ROI_temp','ROI','frqmth','frqs','blnames','measures')

%% Get Time Series

load RSA_DiffFreq

RSA_DiffFreq_Time = [];
RSA_DiffFreq_Time.TimeVec = RSA_DiffFreq.TimeVec;
for r = ROI
    for fm = frqmth
        for fq = frqs
            for bl = blnames
                for m = 1:length(measures)
                    RSA_DiffFreq_Time.(r{1}).(fm{1}).(fq{1}).(bl{1}).(measures{m}).Perceptual.Drawing  = nanmean(RSA_DiffFreq.(r{1}).(fm{1}).(fq{1}).(bl{1}){1,m}(:,:,RSA_DiffFreq.Prcp_Hyp16(:) == 1),3);
                    RSA_DiffFreq_Time.(r{1}).(fm{1}).(fq{1}).(bl{1}).(measures{m}).Perceptual.Picture  = nanmean(RSA_DiffFreq.(r{1}).(fm{1}).(fq{1}).(bl{1}){1,m}(:,:,RSA_DiffFreq.Prcp_Hyp16(:) == 2),3);
                    RSA_DiffFreq_Time.(r{1}).(fm{1}).(fq{1}).(bl{1}).(measures{m}).Perceptual.Within   = nanmean(RSA_DiffFreq.(r{1}).(fm{1}).(fq{1}).(bl{1}){1,m}(:,:,RSA_DiffFreq.Prcp_Hyp16(:) > 0),3);
                    RSA_DiffFreq_Time.(r{1}).(fm{1}).(fq{1}).(bl{1}).(measures{m}).Perceptual.Between  = nanmean(RSA_DiffFreq.(r{1}).(fm{1}).(fq{1}).(bl{1}){1,m}(:,:,RSA_DiffFreq.Prcp_Hyp16(:) < 0),3);
                    RSA_DiffFreq_Time.(r{1}).(fm{1}).(fq{1}).(bl{1}).(measures{m}).Semantic.Animate    = nanmean(RSA_DiffFreq.(r{1}).(fm{1}).(fq{1}).(bl{1}){1,m}(:,:,RSA_DiffFreq.Sem_Hyp16(:) == 1),3);
                    RSA_DiffFreq_Time.(r{1}).(fm{1}).(fq{1}).(bl{1}).(measures{m}).Semantic.NonAnimate = nanmean(RSA_DiffFreq.(r{1}).(fm{1}).(fq{1}).(bl{1}){1,m}(:,:,RSA_DiffFreq.Sem_Hyp16(:) == 2),3);
                    RSA_DiffFreq_Time.(r{1}).(fm{1}).(fq{1}).(bl{1}).(measures{m}).Semantic.Within     = nanmean(RSA_DiffFreq.(r{1}).(fm{1}).(fq{1}).(bl{1}){1,m}(:,:,RSA_DiffFreq.Sem_Hyp16(:) > 0),3);
                    RSA_DiffFreq_Time.(r{1}).(fm{1}).(fq{1}).(bl{1}).(measures{m}).Semantic.Between    = nanmean(RSA_DiffFreq.(r{1}).(fm{1}).(fq{1}).(bl{1}){1,m}(:,:,RSA_DiffFreq.Sem_Hyp16(:) < 0),3);
                end
            end
        end
    end
end




%% Plot RSA Frequency Time series

TimeVec = RSA_DiffFreq_Time.TimeVec;

r = 1;
fm = 3;
fq = 2;
bl = 2;
m = 1;
Cat = {'Perceptual','Semantic'}; c = 1;
dt = 1:2;

figure('Pos', [325 510 650 402]);
plot([]); hold on
Dat_names = fieldnames(RSA_DiffFreq_Time.(ROI{r}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c}));
dat1 = RSA_DiffFreq_Time.(ROI{r}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c}).(Dat_names{dt(1)})(:,:);
dat2 = RSA_DiffFreq_Time.(ROI{r}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c}).(Dat_names{dt(2)})(:,:);
SEM1 = nanstd(dat1(:,:),0,1)./sqrt(size(dat1,1));
SEM2 = nanstd(dat2(:,:),0,1)./sqrt(size(dat2,1));
fill([TimeVec fliplr(TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
h1 = plot(TimeVec, nanmean(dat1,1),'b','linewidth',2);
h2 = plot(TimeVec, nanmean(dat2,1),'r','linewidth',2);
hold off
ylabel(measures{m}); xlabel('Time (s)'); title([Cat{c},': ',Dat_names{dt(1)},' vs ',Dat_names{dt(2)}])
xlim([-0.2 1.5]);lg = legend([h1 h2], {Dat_names{dt(1)},Dat_names{dt(2)}}); legend boxoff; set(lg,'FontSize',14)
box off; %ylim([0.42 0.6])
set(gca,'linewidth',2.5,'FontSize',14)
%saveas(gcf,sprintf('Results/%s_16_%s_LDA_%s.png',ROI{r},Cat{c},datanames{d}))
%close(gcf)



r = 1;
fm = 3;
fq = 2;
bl = 2;
m = 1;
Cat = {'Perceptual','Semantic'}; c = 1;
dt = 3:4;

figure('Pos', [63 66 1802 910]); sbp = 1;
for fq = 1:length(frqs)
    for fm = 1:length(frqmth)
        for bl = 2%:length(blnames)
            subplot(5,3,sbp)
            plot([]); hold on
            Dat_names = fieldnames(RSA_DiffFreq_Time.(ROI{r}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c}));
            dat1 = RSA_DiffFreq_Time.(ROI{r}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c}).(Dat_names{dt(1)})(:,:);
            dat2 = RSA_DiffFreq_Time.(ROI{r}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c}).(Dat_names{dt(2)})(:,:);
            SEM1 = nanstd(dat1(:,:),0,1)./sqrt(size(dat1,1));
            SEM2 = nanstd(dat2(:,:),0,1)./sqrt(size(dat2,1));
            fill([TimeVec fliplr(TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
            fill([TimeVec fliplr(TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
            fill([TimeVec fliplr(TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
            fill([TimeVec fliplr(TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
            h1 = plot(TimeVec, nanmean(dat1,1),'b','linewidth',2);
            h2 = plot(TimeVec, nanmean(dat2,1),'r','linewidth',2);
            hold off
            if(mod(sbp,3) == 1)
                ylabel(frqs{fq});
            end
            if(sbp == 13)
                xlabel('Time (s)');
            end
            if(sbp == 1 || sbp == 2 || sbp == 3)
                title([measures{m},' - ',frqmth{fm},' - ',blnames{bl}])
            end
            if(sbp == 2)
                lg = legend([h1 h2], {[Cat{c},' ',Dat_names{dt(1)}],[Cat{c},' ',Dat_names{dt(2)}]},'location','best'); legend boxoff; set(lg,'FontSize',10)
            end
            xlim([-0.2 1.5]); box off; set(gca,'linewidth',2,'FontSize',10); %ylim([0.42 0.6])
            sbp = sbp + 1;
        end
    end
end


%% Plot RSA Frequency Difference Time series

r1 = 1; r2 = 2;
c1 = 1; c2 = 2;
dt1 = [3 4]; dt2 = [3 4];
Cat = {'Perceptual','Semantic'};
fm = 3;
fq = 5;
bl = 3;
m = 1;

figure('Pos', [325 510 650 402]);
plot([]); hold on
Dat_names1 = fieldnames(RSA_DiffFreq_Time.(ROI{r1}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c1}));
dat1 = RSA_DiffFreq_Time.(ROI{r1}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c1}).(Dat_names1{dt1(2)}) - RSA_DiffFreq_Time.(ROI{r1}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c1}).(Dat_names1{dt1(1)});
Dat_names2 = fieldnames(RSA_DiffFreq_Time.(ROI{r2}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c2}));
dat2 = RSA_DiffFreq_Time.(ROI{r2}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c2}).(Dat_names2{dt2(2)}) - RSA_DiffFreq_Time.(ROI{r2}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c2}).(Dat_names2{dt2(1)});
SEM1 = nanstd(dat1,0,1)./sqrt(size(dat1,1)); SEM2 = nanstd(dat2,0,1)./sqrt(size(dat2,1));
% SEM1 = zeros(1,length(TimeVec)); SEM2 = zeros(1,length(TimeVec));
% for tp = 1:length(TimeVec)
%     [~,~,~,stats] = ttest(RSA_Time.(ROI{r1}).(Cat{c1}).(Dat_names{dt(2)})(:,tp), RSA_Time.(ROI{r1}).(Cat{c1}).(Dat_names{dt(1)})(:,tp));
%     SEM1(1,tp) = tinv(.95,stats.df)*(stats.sd/sqrt(size(dat1,1)));
%     [~,~,~,stats] = ttest(RSA_Time.(ROI{r2}).(Cat{c2}).(Dat_names{dt(2)})(:,tp), RSA_Time.(ROI{r2}).(Cat{c2}).(Dat_names{dt(1)})(:,tp));
%     SEM2(1,tp) = tinv(.95,stats.df)*(stats.sd/sqrt(size(dat2,1)));
% end
fill([TimeVec fliplr(TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
xlim([-0.2 1.5]); %ylim([0.35 0.6]);
plot([-0.2 1.5],[0 0],'--k','linewidth',1)
h1 = plot(TimeVec, nanmean(dat1,1),'b','linewidth',2);
h2 = plot(TimeVec, nanmean(dat2,1),'r','linewidth',2);
ylabel(['Diff ',measures{m}]); xlabel('Time (s)'); %title([Cat{c},': ',Dat_names{dt(1)},' vs ',Dat_names{dt(2)}])
lg = legend([h1 h2], {[ROI{r1},' ',Cat{c1},' ',Dat_names1{dt1(2)},' - ',Dat_names1{dt1(1)}],[ROI{r2},' ',Cat{c2},' ',Dat_names2{dt2(2)},' - ',Dat_names2{dt2(1)}]}); legend boxoff; set(lg,'FontSize',14)
box off;
set(gca,'linewidth',2.5,'FontSize',14,'xlim',[-0.2 1])
hold off
%saveas(gcf,'Results/OCC_TMP_PerceptualvsSemantic_BT-WI_Dat16_LDA.png')
%close(gcf)



r1 = 1; r2 = 1;
c1 = 1; c2 = 2;
dt1 = [3 4]; dt2 = [3 4];
Cat = {'Perceptual','Semantic'};
m = 1;

figure('Pos', [63 66 1802 910]); sbp = 1;
[ha, pos] = tight_subplot(5,3,[.02 .025],[.07 0.05]);
for fq = 1:length(frqs)
    for fm = 1:length(frqmth)
        for bl = 3%:length(blnames)
            axes(ha(sbp)); hold on
            Dat_names1 = fieldnames(RSA_DiffFreq_Time.(ROI{r1}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c1}));
            dat1 = RSA_DiffFreq_Time.(ROI{r1}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c1}).(Dat_names1{dt1(2)}) - RSA_DiffFreq_Time.(ROI{r1}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c1}).(Dat_names1{dt1(1)});
            Dat_names2 = fieldnames(RSA_DiffFreq_Time.(ROI{r2}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c2}));
            dat2 = RSA_DiffFreq_Time.(ROI{r2}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c2}).(Dat_names2{dt2(2)}) - RSA_DiffFreq_Time.(ROI{r2}).(frqmth{fm}).(frqs{fq}).(blnames{bl}).(measures{m}).(Cat{c2}).(Dat_names2{dt2(1)});
            SEM1 = nanstd(dat1,0,1)./sqrt(size(dat1,1)); SEM2 = nanstd(dat2,0,1)./sqrt(size(dat2,1));
            % SEM1 = zeros(1,length(TimeVec)); SEM2 = zeros(1,length(TimeVec));
            % for tp = 1:length(TimeVec)
            %     [~,~,~,stats] = ttest(RSA_Time.(ROI{r1}).(Cat{c1}).(Dat_names{dt(2)})(:,tp), RSA_Time.(ROI{r1}).(Cat{c1}).(Dat_names{dt(1)})(:,tp));
            %     SEM1(1,tp) = tinv(.95,stats.df)*(stats.sd/sqrt(size(dat1,1)));
            %     [~,~,~,stats] = ttest(RSA_Time.(ROI{r2}).(Cat{c2}).(Dat_names{dt(2)})(:,tp), RSA_Time.(ROI{r2}).(Cat{c2}).(Dat_names{dt(1)})(:,tp));
            %     SEM2(1,tp) = tinv(.95,stats.df)*(stats.sd/sqrt(size(dat2,1)));
            % end
            fill([TimeVec fliplr(TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
            fill([TimeVec fliplr(TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
            fill([TimeVec fliplr(TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
            fill([TimeVec fliplr(TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
            xlim([-0.2 1.5]); %ylim([0.35 0.6]);
            plot([-0.2 1.5],[0 0],'--k','linewidth',1)
            h1 = plot(TimeVec, nanmean(dat1,1),'b','linewidth',2);
            h2 = plot(TimeVec, nanmean(dat2,1),'r','linewidth',2);
            set(ha(sbp),'ylim', [-0.06 0.06],'ytick',[-0.03 0 0.03],'yticklabel',{'-0.03' '0' '0.03'})
            if(mod(sbp,3) == 1)
                ylabel(frqs{fq});
            end
            if(sbp == 13)
                xlabel('Time (s)');
                set(ha(sbp),'xticklabel',[-0.2 0 0.2 0.4 0.6 0.8 1])
            end
            if(sbp == 1 || sbp == 2 || sbp == 3)
                title([measures{m},' - ',frqmth{fm},' - ',blnames{bl}])
            end
            if(sbp == 2)
                lg = legend([h1 h2], {[ROI{r1},' ',Cat{c1},' ',Dat_names1{dt1(2)},' - ',Dat_names1{dt1(1)}],[ROI{r2},' ',Cat{c2},' ',Dat_names2{dt2(2)},' - ',Dat_names2{dt2(1)}]}); 
                legend boxoff; set(lg,'FontSize',10,'location','best')
            end
            box off; set(ha(sbp),'linewidth',2.5,'FontSize',14,'xlim',[-0.2 1])
            hold off
            sbp = sbp + 1;
        end
    end
end


