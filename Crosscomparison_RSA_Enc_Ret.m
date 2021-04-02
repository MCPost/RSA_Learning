%% Cross-comparison of RSA time courses

% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs

% Hypotheses Matrix
trl_mat = [kron([1;2],ones(64,1)) kron([1;2;1;2],ones(32,1))];
Perceptual_Mat_full = zeros(size(trl_mat,1));
Semantic_Mat_full = zeros(size(trl_mat,1));
for i = 1:size(trl_mat,1)-1
    for j = (i+1):size(trl_mat,1)
        if(j ~= size(trl_mat,1) - (i - 1))
            if(trl_mat(i,1) == 1 && trl_mat(j,1) == 1)
                Perceptual_Mat_full(i,j) = 1;
            elseif(trl_mat(i,1) == 2 && trl_mat(j,1) == 2)
                Perceptual_Mat_full(i,j) = 2;
            else
                Perceptual_Mat_full(i,j) = -1;
            end

            if(trl_mat(i,2) == 1 && trl_mat(j,2) == 1)
                Semantic_Mat_full(i,j) = 1;
            elseif(trl_mat(i,2) == 2 && trl_mat(j,2) == 2)
                Semantic_Mat_full(i,j) = 2;
            else
                Semantic_Mat_full(i,j) = -1;
            end
        end
    end
end

Perceptual_Mat_red16 = zeros(16);
Semantic_Mat_red16 = zeros(16);
for i = 1:16-1
    for j = (i+1):16
        if(j ~= 16 - (i - 1))
            if(trl_mat(8*i,1) == 1 && trl_mat(8*j,1) == 1)
                Perceptual_Mat_red16(i,j) = 1;
            elseif(trl_mat(8*i,1) == 2 && trl_mat(8*j,1) == 2)
                Perceptual_Mat_red16(i,j) = 2;
            else
                Perceptual_Mat_red16(i,j) = -1;
            end

            if(trl_mat(8*i,2) == 1 && trl_mat(8*j,2) == 1)
                Semantic_Mat_red16(i,j) = 1;
            elseif(trl_mat(8*i,2) == 2 && trl_mat(8*j,2) == 2)
                Semantic_Mat_red16(i,j) = 2;
            else
                Semantic_Mat_red16(i,j) = -1;
            end
        end
    end
end



%% Load Data 

% Subject Names
load('RSA_Data_Enc','Subj_names')

% Measure
msr = 1;

%% Cross-correlate RSA time courses for each subject

% save('CrossComp_RSA','Subj_names')
% 
% c = parcluster();
% jobHandles = cell(4,1);
% for msr = [1 3]
%     jobHandles{msr} = batch(c, @Crosscomparison_RSA_Enc_Ret_LoadData, 0, {msr});
% end
% 
% jobHandles{1}
% diary(jobHandles{1})

% cfg = [];
% cfg.slide_window   = 0.050;
% cfg.slide_step     = 0.010;
% cfg.window_average = 'gaussian';
% cfg.Hyp_perceptual = {Perceptual_Mat_full Perceptual_Mat_red16};
% cfg.Hyp_semantic   = {Semantic_Mat_full   Semantic_Mat_red16};
% cfg.ROI            = {'OCC','TMP','FRT','CNT','PRT'};
% cfg.only16         = true;
% CrossComp_RSA = create_xcomp_rsa(cfg, RSA_Data_Enc, RSA_Data_Ret);


load('CrossComp_RSA')


Data = CrossComp_RSA_LDA.RSA_red16;
zmapthresh1 = CrossComp_RSA_LDA.zmapthresh{1}; zmapthresh1(isnan(zmapthresh1)) = 0;zmapthresh1(zmapthresh1 ~= 0) = 1;
zmapthresh2 = CrossComp_RSA_LDA.zmapthresh{2}; zmapthresh2(isnan(zmapthresh2)) = 0;zmapthresh2(zmapthresh2 ~= 0) = 1;
%zmapthresh3 = CrossComp_RSA_LDA.zmapthresh{3}; zmapthresh3(isnan(zmapthresh3)) = 0;zmapthresh3(zmapthresh3 ~= 0) = 1;
%zmapthresh4 = CrossComp_RSA_LDA.zmapthresh{4}; zmapthresh4(isnan(zmapthresh4)) = 0;zmapthresh4(zmapthresh4 ~= 0) = 1;
%zmapthresh5 = CrossComp_RSA_LDA.zmapthresh{5}; zmapthresh5(isnan(zmapthresh5)) = 0;zmapthresh5(zmapthresh5 ~= 0) = 1;
TimeX = CrossComp_RSA_LDA.TimeVec1;
TimeY = CrossComp_RSA_LDA.TimeVec2;


figure
subplot(1,2,1)
contourf(TimeX, TimeY, squeeze(nanmean(Data.OCC,1)), 40,'linestyle','none'); colorbar
caxis([-.07 .07]); set(gca,'xlim', [-0.2 1.5], 'ylim', [-4 0.2]); title('OCCIPITAL'); xlabel('Encoding'); ylabel('Retrieval')
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
contour(TimeX, TimeY,zmapthresh1,1,'linecolor','k','linewidth',1.5)
hold off
subplot(1,2,2)
contourf(TimeX, TimeY, squeeze(nanmean(Data.TMP,1)), 40,'linestyle','none'); colorbar
caxis([-.07 .07]); set(gca,'xlim', [-0.2 1.5], 'ylim', [-4 0.2]); title('TEMPORAL'); xlabel('Encoding'); ylabel('Retrieval')
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
contour(TimeX, TimeY,zmapthresh2,1,'linecolor','k','linewidth',1.5)
hold off



TimeX = CrossComp_RSA_euclidian_wcc.TimeVec1;
TimeY = CrossComp_RSA_euclidian_wcc.TimeVec2;
enc_win = [-0.2 0.6];
ret_win = [-2.5 0.2];

figure
subplot(2,4,1)
contourf(TimeX, TimeY, squeeze(nanmean(CrossComp_RSA_LDA.RSA_red16.OCC,1)), 40,'linestyle','none'); colorbar
caxis([-.07 .07]); set(gca,'xlim', enc_win, 'ylim', ret_win); title('OCCIPITAL'); xlabel('Encoding'); ylabel('Retrieval')
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
zmapthresh = CrossComp_RSA_LDA.zmapthresh{1}; zmapthresh(isnan(zmapthresh)) = 0;zmapthresh(zmapthresh ~= 0) = 1;
contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
hold off
subplot(2,4,5)
contourf(TimeX, TimeY, squeeze(nanmean(CrossComp_RSA_LDA.RSA_red16.TMP,1)), 40,'linestyle','none'); colorbar
caxis([-.07 .07]); set(gca,'xlim', enc_win, 'ylim', ret_win); title('TEMPORAL'); xlabel('Encoding'); ylabel('Retrieval')
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
zmapthresh = CrossComp_RSA_LDA.zmapthresh{2}; zmapthresh(isnan(zmapthresh)) = 0;zmapthresh(zmapthresh ~= 0) = 1;
contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
hold off

subplot(2,4,2)
contourf(TimeX, TimeY, squeeze(nanmean(CrossComp_RSA_SVM.RSA_red16.OCC,1)), 40,'linestyle','none'); colorbar
caxis([-.07 .07]); set(gca,'xlim', enc_win, 'ylim', ret_win); title('OCCIPITAL'); xlabel('Encoding'); ylabel('Retrieval')
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
zmapthresh = CrossComp_RSA_SVM.zmapthresh{1}; zmapthresh(isnan(zmapthresh)) = 0;zmapthresh(zmapthresh ~= 0) = 1;
contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
hold off
subplot(2,4,6)
contourf(TimeX, TimeY, squeeze(nanmean(CrossComp_RSA_SVM.RSA_red16.TMP,1)), 40,'linestyle','none'); colorbar
caxis([-.07 .07]); set(gca,'xlim', enc_win, 'ylim', ret_win); title('TEMPORAL'); xlabel('Encoding'); ylabel('Retrieval')
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
zmapthresh = CrossComp_RSA_SVM.zmapthresh{2}; zmapthresh(isnan(zmapthresh)) = 0;zmapthresh(zmapthresh ~= 0) = 1;
contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
hold off

subplot(2,4,3)
contourf(TimeX, TimeY, squeeze(nanmean(CrossComp_RSA_euclidian.RSA_red16.OCC,1)), 40,'linestyle','none'); colorbar
caxis([-.07 .07]); set(gca,'xlim', enc_win, 'ylim', ret_win); title('OCCIPITAL'); xlabel('Encoding'); ylabel('Retrieval')
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
zmapthresh = CrossComp_RSA_euclidian.zmapthresh{1}; zmapthresh(isnan(zmapthresh)) = 0;zmapthresh(zmapthresh ~= 0) = 1;
contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
hold off
subplot(2,4,7)
contourf(TimeX, TimeY, squeeze(nanmean(CrossComp_RSA_euclidian.RSA_red16.TMP,1)), 40,'linestyle','none'); colorbar
caxis([-.07 .07]); set(gca,'xlim', enc_win, 'ylim', ret_win); title('TEMPORAL'); xlabel('Encoding'); ylabel('Retrieval')
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
zmapthresh = CrossComp_RSA_euclidian.zmapthresh{2}; zmapthresh(isnan(zmapthresh)) = 0;zmapthresh(zmapthresh ~= 0) = 1;
contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
hold off

subplot(2,4,4)
contourf(TimeX, TimeY, squeeze(nanmean(CrossComp_RSA_euclidian_wcc.RSA_red16.OCC,1)), 40,'linestyle','none'); colorbar
caxis([-.07 .07]); set(gca,'xlim', enc_win, 'ylim', ret_win); title('OCCIPITAL'); xlabel('Encoding'); ylabel('Retrieval')
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
zmapthresh = CrossComp_RSA_euclidian_wcc.zmapthresh{1}; zmapthresh(isnan(zmapthresh)) = 0;zmapthresh(zmapthresh ~= 0) = 1;
contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
hold off
subplot(2,4,8)
contourf(TimeX, TimeY, squeeze(nanmean(CrossComp_RSA_euclidian_wcc.RSA_red16.TMP,1)), 40,'linestyle','none'); colorbar
caxis([-.07 .07]); set(gca,'xlim', enc_win, 'ylim', ret_win); title('TEMPORAL'); xlabel('Encoding'); ylabel('Retrieval')
hold on
plot([0 0],[TimeY(1) TimeY(end)],'--w','linewidth',2)
plot([TimeX(1) TimeX(end)],[0 0],'--w','linewidth',2)
zmapthresh = CrossComp_RSA_euclidian_wcc.zmapthresh{2}; zmapthresh(isnan(zmapthresh)) = 0;zmapthresh(zmapthresh ~= 0) = 1;
contour(TimeX, TimeY,zmapthresh,1,'linecolor','k','linewidth',1.5)
hold off



