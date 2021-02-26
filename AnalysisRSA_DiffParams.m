%% Different Parameter Permutations Analysis

%% Subject Names
Subj_names = {'AG','CEWD','CM','DS','FVM'};

% %% Get Files
% mat_path = 'RSA_DiffParams/';
% load([mat_path,'Elec_Info'])
% 
% % Get Filenames
% filenames = dir(mat_path);
% filenames = struct2cell(filenames);
% filenames = filenames(1,:)';
% filenames = filenames(~cellfun(@isempty, regexp(filenames,'RSA')) & ~cellfun(@isempty, regexp(filenames,'Data_CM')));
% 
% % Create Data Struct
% RSA_DiffParams = []; datanames = {};
% for i = 1:length(filenames)
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat128.Perceptual.Drawing = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat128.Perceptual.Picture = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat128.Perceptual.Within = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat128.Perceptual.Between = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat128.Semantic.Animate = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat128.Semantic.NonAnimate = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat128.Semantic.Within = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat128.Semantic.Between = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat16.Perceptual.Drawing = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat16.Perceptual.Picture = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat16.Perceptual.Within = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat16.Perceptual.Between = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat16.Semantic.Animate = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat16.Semantic.NonAnimate = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat16.Semantic.Within = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).OCC.Dat16.Semantic.Between = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat128.Perceptual.Drawing = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat128.Perceptual.Picture = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat128.Perceptual.Within = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat128.Perceptual.Between = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat128.Semantic.Animate = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat128.Semantic.NonAnimte = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat128.Semantic.Within = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat128.Semantic.Between = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat16.Perceptual.Drawing = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat16.Perceptual.Picture = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat16.Perceptual.Within = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat16.Perceptual.Between = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat16.Semantic.Animate = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat16.Semantic.NonAnimate = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat16.Semantic.Within = nan(5,5,436);
%     RSA_DiffParams.(filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4)).TMP.Dat16.Semantic.Between = nan(5,5,436);
%     datanames{i,1} = filenames{i}(regexp(filenames{1}, '_ArtCorr')+1:end-4);
% end
% 
% % Hypotheses Matrix
% trl_mat = [kron([1;2],ones(64,1)) kron([1;2;1;2],ones(32,1))];
% RSA_DiffParams.Prcp_Hyp128 = zeros(size(trl_mat,1));
% RSA_DiffParams.Sem_Hyp128 = zeros(size(trl_mat,1));
% for i = 1:size(trl_mat,1)-1
%     for j = (i+1):size(trl_mat,1)
%         if(j ~= size(trl_mat,1) - (i - 1))
%             if(trl_mat(i,1) == 1 && trl_mat(j,1) == 1)
%                 RSA_DiffParams.Prcp_Hyp128(i,j) = 1;
%             elseif(trl_mat(i,1) == 2 && trl_mat(j,1) == 2)
%                 RSA_DiffParams.Prcp_Hyp128(i,j) = 2;
%             else
%                 RSA_DiffParams.Prcp_Hyp128(i,j) = -1;
%             end
% 
%             if(trl_mat(i,2) == 1 && trl_mat(j,2) == 1)
%                 RSA_DiffParams.Sem_Hyp128(i,j) = 1;
%             elseif(trl_mat(i,2) == 2 && trl_mat(j,2) == 2)
%                 RSA_DiffParams.Sem_Hyp128(i,j) = 2;
%             else
%                 RSA_DiffParams.Sem_Hyp128(i,j) = -1;
%             end
%         end
%     end
% end
% 
% RSA_DiffParams.Prcp_Hyp16 = zeros(16);
% RSA_DiffParams.Sem_Hyp16 = zeros(16);
% for i = 1:16-1
%     for j = (i+1):16
%         if(j ~= 16 - (i - 1))
%             if(trl_mat(8*i,1) == 1 && trl_mat(8*j,1) == 1)
%                 RSA_DiffParams.Prcp_Hyp16(i,j) = 1;
%             elseif(trl_mat(8*i,1) == 2 && trl_mat(8*j,1) == 2)
%                 RSA_DiffParams.Prcp_Hyp16(i,j) = 2;
%             else
%                 RSA_DiffParams.Prcp_Hyp16(i,j) = -1;
%             end
% 
%             if(trl_mat(8*i,2) == 1 && trl_mat(8*j,2) == 1)
%                 RSA_DiffParams.Sem_Hyp16(i,j) = 1;
%             elseif(trl_mat(8*i,2) == 2 && trl_mat(8*j,2) == 2)
%                 RSA_DiffParams.Sem_Hyp16(i,j) = 2;
%             else
%                 RSA_DiffParams.Sem_Hyp16(i,j) = -1;
%             end
%         end
%     end
% end
% clear trl_mat
% 
% 
% %% Import Data
% 
% for sub = 1:length(Subj_names)
%     
%     fprintf('\n')
%     fprintf('Subject %s (%i / %i) started loading!', Subj_names{sub}, sub, length(Subj_names))
%     fprintf('\n')
%     
%     fprintf('\n')
%     fprintf('Subject: %s  ---  Dataset Loading Progress:  0.0 %%',Subj_names{sub})
%     txt = sprintf('Subject: %s  ---  Dataset Loading Progress:  0.0 %%',Subj_names{sub});
%     for d = 1:length(datanames)
%         temp = load([mat_path,'RSA_Data_',Subj_names{sub},'_',datanames{d}]);
%         temp_names = fieldnames(temp);
%         if(length(temp_names) > 2)
%             OCC = temp.(temp_names{2}); TMP = temp.(temp_names{3});
%             clear temp temp_names
% 
%             for tp = 1:length(OCC.TimeVec)
% 
%                 if(isempty(regexp(datanames{d},'NoMNN','ONCE')))
% 
%                     curd = atanh(OCC.RSA_full{1,1}(:,:,tp));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Drawing(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Picture(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Within(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Between(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 < 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Animate(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.NonAnimate(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Within(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Between(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 < 0));
%                     curd = atanh(OCC.RSA_full{1,2}(:,:,tp));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Drawing(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Picture(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Within(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Between(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 < 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Animate(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.NonAnimate(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Within(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Between(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 < 0));
%                     curd = OCC.RSA_full{1,3}(:,:,tp);
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Drawing(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Picture(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Within(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Between(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 < 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Animate(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.NonAnimate(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Within(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Between(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 < 0));
%                     curd = OCC.RSA_full{1,4}(:,:,tp);
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Drawing(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Picture(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Within(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Between(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 < 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Animate(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.NonAnimate(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Within(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Between(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 < 0));
%                     curd = atanh(OCC.RSA_full{1,5}(:,:,tp));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Drawing(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Picture(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Within(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Perceptual.Between(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 < 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Animate(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.NonAnimate(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Within(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).OCC.Dat128.Semantic.Between(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 < 0));
% 
%                     curd = OCC.RSA_full{1,1}(:,:,tp);
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Drawing(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Picture(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Within(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Between(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 < 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Animate(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.NonAnimate(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Within(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Between(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 < 0));
%                     curd = TMP.RSA_full{1,2}(:,:,tp);
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Drawing(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Picture(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Within(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Between(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 < 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Animate(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.NonAnimate(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Within(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Between(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 < 0));
%                     curd = TMP.RSA_full{1,3}(:,:,tp);
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Drawing(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Picture(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Within(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Between(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 < 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Animate(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.NonAnimate(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Within(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Between(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 < 0));
%                     curd = TMP.RSA_full{1,4}(:,:,tp);
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Drawing(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Picture(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Within(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Between(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 < 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Animate(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.NonAnimate(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Within(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Between(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 < 0));
%                     curd = TMP.RSA_full{1,5}(:,:,tp);
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Drawing(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Picture(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Within(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Perceptual.Between(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp128 < 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Animate(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 1));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.NonAnimate(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 == 2));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Within(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 > 0));
%                     RSA_DiffParams.(datanames{d}).TMP.Dat128.Semantic.Between(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp128 < 0));
% 
%                 end
% 
%                 curd = OCC.RSA_16{1,1}(:,:,tp);
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Drawing(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Picture(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Within(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Between(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 < 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Animate(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.NonAnimate(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Within(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Between(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 < 0));
%                 curd = OCC.RSA_16{1,2}(:,:,tp);
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Drawing(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Picture(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Within(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Between(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 < 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Animate(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.NonAnimate(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Within(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Between(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 < 0));
%                 curd = OCC.RSA_16{1,3}(:,:,tp);
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Drawing(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Picture(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Within(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Between(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 < 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Animate(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.NonAnimate(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Within(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Between(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 < 0));
%                 curd = OCC.RSA_16{1,4}(:,:,tp);
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Drawing(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Picture(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Within(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Between(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 < 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Animate(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.NonAnimate(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Within(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Between(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 < 0));
%                 curd = OCC.RSA_16{1,5}(:,:,tp);
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Drawing(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Picture(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Within(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Perceptual.Between(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 < 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Animate(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.NonAnimate(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Within(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).OCC.Dat16.Semantic.Between(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 < 0));
% 
%                 curd = OCC.RSA_16{1,1}(:,:,tp);
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Drawing(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Picture(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Within(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Between(1,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 < 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Animate(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.NonAnimate(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Within(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Between(1,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 < 0));
%                 curd = TMP.RSA_16{1,2}(:,:,tp);
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Drawing(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Picture(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Within(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Between(2,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 < 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Animate(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.NonAnimate(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Within(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Between(2,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 < 0));
%                 curd = TMP.RSA_16{1,3}(:,:,tp);
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Drawing(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Picture(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Within(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Between(3,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 < 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Animate(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.NonAnimate(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Within(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Between(3,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 < 0));
%                 curd = TMP.RSA_16{1,4}(:,:,tp);
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Drawing(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Picture(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Within(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Between(4,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 < 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Animate(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.NonAnimate(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Within(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Between(4,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 < 0));
%                 curd = TMP.RSA_16{1,5}(:,:,tp);
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Drawing(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Picture(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Within(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Perceptual.Between(5,sub,tp) = nanmean(curd(RSA_DiffParams.Prcp_Hyp16 < 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Animate(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 1));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.NonAnimate(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 == 2));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Within(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 > 0));
%                 RSA_DiffParams.(datanames{d}).TMP.Dat16.Semantic.Between(5,sub,tp) = nanmean(curd(RSA_DiffParams.Sem_Hyp16 < 0));
% 
% 
%             end
%             
%         end
%         
%         fprintf(repmat('\b',1,length(txt)))
%         fprintf('Subject: %s  ---  Dataset Loading Progress: %3.2f %%',Subj_names{sub},(d/length(datanames))*100)
%         txt = sprintf('Subject: %s  ---  Dataset Loading Progress: %3.2f %%',Subj_names{sub},(d/length(datanames))*100);
%     end
%     
%     fprintf('\n')
%     fprintf('Subject %s (%i / %i) finished loading!', Subj_names{sub}, sub, length(Subj_names))
%     fprintf('\n')
%     
% end
% RSA_DiffParams.TimeVec = OCC.TimeVec;
% RSA_DiffParams.Meas128 = {'pearson' 'spearman' 'euclidian' 'mahalanobis' 'cosine'};
% RSA_DiffParams.Meas16 = {'euclidian c.v.' 'pearson c.v.' 'LDA' 'SVM' 'DistCorr'};
% 
% save('RSA_DiffParams','RSA_DiffParams','Elec_Info','Subj_names')
% 


%% Plot Results

load RSA_DiffParams

datanames = fieldnames(RSA_DiffParams);
datanames = datanames(1:24);

ROI = {'OCC','TMP'}; r = 1;
MSize = {'Dat128','Dat16'}; sz = 1;
Cat = {'Perceptual','Semantic'}; c = 1;
dt = 3:4;
Measure = {'Meas128','Meas16'}; m = 2;

for d = [5 6 7 8] %1:length(datanames)
    for r = 2%:length(ROI)
        for sz = 2%:length(MSize)
            for c = 2%:length(Cat)
                for m = 4
                    figure('Pos', [325 510 650 402]);
                    plot([]); hold on
                    Dat_names = fieldnames(RSA_DiffParams.(datanames{d}).(ROI{r}).(MSize{sz}).(Cat{c}));
                    dat1 = squeeze(RSA_DiffParams.(datanames{d}).(ROI{r}).(MSize{sz}).(Cat{c}).(Dat_names{dt(1)})(m,:,:));
                    dat2 = squeeze(RSA_DiffParams.(datanames{d}).(ROI{r}).(MSize{sz}).(Cat{c}).(Dat_names{dt(2)})(m,:,:));
                    fill([RSA_DiffParams.TimeVec fliplr(RSA_DiffParams.TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) + nanstd(dat1,1)./sqrt(size(dat1,1)))],'b','FaceAlpha',0.3,'EdgeAlpha',0);
                    fill([RSA_DiffParams.TimeVec fliplr(RSA_DiffParams.TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) - nanstd(dat1,1)./sqrt(size(dat1,1)))],'b','FaceAlpha',0.3,'EdgeAlpha',0);
                    fill([RSA_DiffParams.TimeVec fliplr(RSA_DiffParams.TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) + nanstd(dat2,1)./sqrt(size(dat2,1)))],'r','FaceAlpha',0.3,'EdgeAlpha',0);
                    fill([RSA_DiffParams.TimeVec fliplr(RSA_DiffParams.TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) - nanstd(dat2,1)./sqrt(size(dat2,1)))],'r','FaceAlpha',0.3,'EdgeAlpha',0);
                    h1 = plot(RSA_DiffParams.TimeVec, nanmean(dat1,1),'b','linewidth',2);
                    h2 = plot(RSA_DiffParams.TimeVec, nanmean(dat2,1),'r','linewidth',2);
                    hold off
                    ylabel(RSA_DiffParams.(Measure{sz}){m}); xlabel('Time (s)'); title([Cat{c},': ',Dat_names{dt(1)},' vs ',Dat_names{dt(2)}])
                    xlim([-0.2 1.5]);lg = legend([h1 h2], {Dat_names{dt(1)},Dat_names{dt(2)}}); legend boxoff; set(lg,'FontSize',14)
                    box off; %ylim([0.35 0.6])
                    set(gca,'linewidth',2.5,'FontSize',14)
                    %saveas(gcf,sprintf('Results/DiffParams/%s_%s_%s_%s_%s.png',ROI{r},MSize{sz},Cat{c},RSA_DiffParams.(Measure{sz}){m},datanames{d}))
                    %close(gcf)
                end
            end
        end
    end
end



r1 = 1; r2 = 2;
c1 = 1; c2 = 2;
dt = [3 4];
for d = [6] %1:length(datanames)
    for sz = 2%:length(MSize)
        for m = 3
            figure('Pos', [325 510 650 402]);
            plot([]); hold on
            Dat_names = fieldnames(RSA_DiffParams.(datanames{d}).(ROI{r}).(MSize{sz}).(Cat{c}));
            dat1 = squeeze(RSA_DiffParams.(datanames{d}).(ROI{r1}).(MSize{sz}).(Cat{c1}).(Dat_names{dt(2)})(m,:,:)) - squeeze(RSA_DiffParams.(datanames{d}).(ROI{r1}).(MSize{sz}).(Cat{c1}).(Dat_names{dt(1)})(m,:,:));
            dat2 = squeeze(RSA_DiffParams.(datanames{d}).(ROI{r2}).(MSize{sz}).(Cat{c2}).(Dat_names{dt(2)})(m,:,:)) - squeeze(RSA_DiffParams.(datanames{d}).(ROI{r2}).(MSize{sz}).(Cat{c2}).(Dat_names{dt(1)})(m,:,:));
            SEM1 = nanstd(dat1,1)./sqrt(size(dat1,1)); SEM2 = nanstd(dat2,1)./sqrt(size(dat2,1));
%             SEM1 = zeros(1,length(RSA_DiffParams.TimeVec)); SEM2 = zeros(1,length(RSA_DiffParams.TimeVec));
%             for tp = 1:length(RSA_DiffParams.TimeVec)
%                 [~,~,~,stats] = ttest(squeeze(RSA_DiffParams.(datanames{d}).(ROI{r1}).(MSize{sz}).(Cat{c1}).(Dat_names{dt(2)})(m,:,tp))', squeeze(RSA_DiffParams.(datanames{d}).(ROI{r1}).(MSize{sz}).(Cat{c1}).(Dat_names{dt(1)})(m,:,tp))');
%                 SEM1(1,tp) = tinv(.95,stats.df)*(stats.sd/sqrt(size(dat1,1)));
%                 [~,~,~,stats] = ttest(squeeze(RSA_DiffParams.(datanames{d}).(ROI{r2}).(MSize{sz}).(Cat{c2}).(Dat_names{dt(2)})(m,:,tp))', squeeze(RSA_DiffParams.(datanames{d}).(ROI{r2}).(MSize{sz}).(Cat{c2}).(Dat_names{dt(1)})(m,:,tp))');
%                 SEM2(1,tp) = tinv(.95,stats.df)*(stats.sd/sqrt(size(dat2,1)));
%             end
            fill([RSA_DiffParams.TimeVec fliplr(RSA_DiffParams.TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
            fill([RSA_DiffParams.TimeVec fliplr(RSA_DiffParams.TimeVec)],[nanmean(dat1,1) fliplr(nanmean(dat1,1) - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
            fill([RSA_DiffParams.TimeVec fliplr(RSA_DiffParams.TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
            fill([RSA_DiffParams.TimeVec fliplr(RSA_DiffParams.TimeVec)],[nanmean(dat2,1) fliplr(nanmean(dat2,1) - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
            xlim([-0.1 1.5]); %ylim([0.35 0.6]);
            plot([-0.1 1.5],[0 0],'--k','linewidth',1)
            h1 = plot(RSA_DiffParams.TimeVec, nanmean(dat1,1),'b','linewidth',2);
            h2 = plot(RSA_DiffParams.TimeVec, nanmean(dat2,1),'r','linewidth',2);
            hold off
            ylabel(RSA_DiffParams.(Measure{sz}){m}); xlabel('Time (s)'); %title([Cat{c},': ',Dat_names{dt(1)},' vs ',Dat_names{dt(2)}])
            lg = legend([h1 h2], {[ROI{r1},' ',Cat{c1},' BT - WI'],[ROI{r2},' ',Cat{c2},' BT - WI']}); legend boxoff; set(lg,'FontSize',14)
            box off;
            set(gca,'linewidth',2.5,'FontSize',14)
            %saveas(gcf,sprintf('Results/DiffParams/%s_%s_%s_%s_%s.png',ROI{r},MSize{sz},Cat{c},RSA_DiffParams.(Measure{sz}){m},datanames{d}))
            %close(gcf)
            saveas(gcf,'Results/DiffParams/OCC_TMP_PerceptualvsSemantic_BT-WI_Dat16_LDA.png')
            close(gcf)
        end
    end
end
