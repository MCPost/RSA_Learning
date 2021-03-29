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


figure
subplot(1,2,1)
contourf(RSA_Data_Enc.TimeVec, RSA_Data_Ret.TimeVec, flip(squeeze(nanmean(CrossComp_RSA.OCC(:,2,:,:) - CrossComp_RSA.OCC(:,1,:,:),1)),1), 40,'linestyle','none')
caxis([-.06 .06]); set(gca,'xlim', [-0.2 1], 'ylim', [-4 0.2])
hold on
plot([0 0],[RSA_Data_Ret.TimeVec(1) RSA_Data_Ret.TimeVec(end)],'--w','linewidth',2)
plot([RSA_Data_Enc.TimeVec(1) RSA_Data_Enc.TimeVec(end)],[0 0],'--w','linewidth',2)
hold off
subplot(1,2,2)
contourf(RSA_Data_Enc.TimeVec, RSA_Data_Ret.TimeVec, flip(squeeze(nanmean(CrossComp_RSA.TMP(:,4,:,:) - CrossComp_RSA.TMP(:,3,:,:),1)),1), 40,'linestyle','none')
caxis([-.06 .06]); set(gca,'xlim', [-0.2 1], 'ylim', [-4 0.2])
hold on
plot([0 0],[RSA_Data_Ret.TimeVec(1) RSA_Data_Ret.TimeVec(end)],'--w','linewidth',2)
plot([RSA_Data_Enc.TimeVec(1) RSA_Data_Enc.TimeVec(end)],[0 0],'--w','linewidth',2)
hold off



