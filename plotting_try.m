%% Plotting Try

Name = 'AG';

% Biosemi 128 Electrodes System Radial ABC
elecs = {'A1';'A2';'A3';'A4';'A5';'A6';'A7';'A8';'A9';'A10';'A11';'A12';'A13';'A14';'A15';'A16';'A17';...
          'A18';'A19';'A20';'A21';'A22';'A23';'A24';'A25';'A26';'A27';'A28';'A29';'A30';'A31';'A32';'B1';...
          'B2';'B3';'B4';'B5';'B6';'B7';'B8';'B9';'B10';'B11';'B12';'B13';'B14';'B15';'B16';'B17';'B18';...
          'B19';'B20';'B21';'B22';'B23';'B24';'B25';'B26';'B27';'B28';'B29';'B30';'B31';'B32';'C1';'C2';...
          'C3';'C4';'C5';'C6';'C7';'C8';'C9';'C10';'C11';'C12';'C13';'C14';'C15';'C16';'C17';'C18';'C19';...
          'C20';'C21';'C22';'C23';'C24';'C25';'C26';'C27';'C28';'C29';'C30';'C31';'C32';'D1';'D2';'D3';...
          'D4';'D5';'D6';'D7';'D8';'D9';'D10';'D11';'D12';'D13';'D14';'D15';'D16';'D17';'D18';'D19';'D20';...
          'D21';'D22';'D23';'D24';'D25';'D26';'D27';'D28';'D29';'D30';'D31';'D32'};
ROI_occ = {'A8','A9','A10','A11','A12','A13','A14','A15','A16','A17','A21','A22','A23','A24','A25','A26',...
           'A27','A28','A29','A30','B5','B6','B7','B8','B9'};
ROI_occ_idx = find(cell2mat(cellfun(@(x) any(strcmp(x, ROI_occ)), elecs, 'UniformOutput', 0)));
curROI = ROI_occ_idx;

Preproc_Data_1 = load(['Prepro_2_encoding_full_',Name,'.mat'], ['r_1_',Name]);
Preproc_Data_1 = Preproc_Data_1.(cell2mat(fieldnames(Preproc_Data_1)));
Preproc_Data_2 = load(['Prepro_2_encoding_full_',Name,'.mat'], ['r_2_',Name]);
Preproc_Data_2 = Preproc_Data_2.(cell2mat(fieldnames(Preproc_Data_2)));

logfile = importdata(['Logs/resultfile_',Name,'.csv']);

logfile.textdata(2:end,find(strcmp(logfile.textdata(find(strcmp(logfile.textdata(:,2),'retrieval'),1,'first'),:),''),1,'first')+(0:size(logfile.data,2)-1)) = num2cell(logfile.data);
logfile = logfile.textdata;
for i = 1:size(logfile,2)
    if(~isempty(str2double(logfile{2,i})) && ~isnan(str2double(logfile{2,i})))
        logfile(2:end,i) = num2cell(str2double(logfile(2:end,i)));
    end
end
clear i


%% Get Encoding from logfile and Categories of Objects

names = {'block_number','trial_encoding','trial_total_enc','object_name','object_number','object_percep',... 
         'object_category','RT_encoding','cue_retrieval','OnsetEncodingBlockStart','WordOnsetEncoding',... 
         'ObjectContextOnset','ResponseEncodingOnset'};

Encoding_Data = logfile(strcmp(logfile(:,2),'encoding'), cell2mat(cellfun(@(x) ismember(x,names), strtrim(logfile(1,:)),'UniformOutput', 0)));
% block_number, trial_encoding, trial_total_enc, object_name,
% object_number, object_percep, object_category, RT_encoding,
% cue_retrieval, OnsetEncodingBlockStart, WordOnsetEncoding, 
% ObjectContextOnset, ResponseEncodingOnset, PerceptualDim, SemanticDim,
% Category

Encoding_Data(cell2mat(Encoding_Data(:,6)) == 1,14) = {'Drawing'};
Encoding_Data(cell2mat(Encoding_Data(:,6)) == 2,14) = {'Picture'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) < 5,15) = {'Animate'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) > 4,15) = {'Inanimate'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) < 5,16) = {1};
Encoding_Data(cell2mat(Encoding_Data(:,7)) > 4,16) = {2};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 1,17) = {'insect'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 2,17) = {'bird'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 3,17) = {'mammal'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 4,17) = {'sea'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 5,17) = {'fruit'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 6,17) = {'electronic'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 7,17) = {'veggie'};
Encoding_Data(cell2mat(Encoding_Data(:,7)) == 8,17) = {'clothe'};


%% Marry Behavioral and EEG Data

TrialInfo = cell(size(Encoding_Data,1),12);
Preproc_Data_2.trialinfo = num2cell(Preproc_Data_2.trialinfo);
draw_idx = find(cell2mat(Encoding_Data(:,6)) == 1); cnt = 1;
for i = 1:length(draw_idx)
    if(ismember(draw_idx(i), cell2mat(Preproc_Data_2.trialinfo(:,2))))
        TrialInfo(i,:) = [Preproc_Data_2.trialinfo(cnt,1:4) Encoding_Data(draw_idx(i),[6 14 16 15 7 17 5 4])];
        cnt = cnt + 1;
    else
        TrialInfo(i,:) = [{0}, {draw_idx(i)} Encoding_Data(draw_idx(i),[6 7 6 14 16 15 7 17 5 4])];
    end
end
Preproc_Data_1.trialinfo = num2cell(Preproc_Data_1.trialinfo);
pic_idx = find(cell2mat(Encoding_Data(:,6)) == 2); cnt = 1;
for i = 1:length(pic_idx)
    if(ismember(pic_idx(i), cell2mat(Preproc_Data_1.trialinfo(:,2))))
        TrialInfo(64+i,:) = [Preproc_Data_1.trialinfo(cnt,1:4) Encoding_Data(pic_idx(i),[6 14 16 15 7 17 5 4])];
        cnt = cnt + 1;
    else
        TrialInfo(64+i,:) = [{0}, {pic_idx(i)} Encoding_Data(pic_idx(i),[6 7 6 14 16 15 7 17 5 4])];
    end
end

% Get Index of NaN trials
Nan_trial_idx = find(cell2mat(TrialInfo(:,1)) == 0);

%% Create EEG Data Matrix
[~,sort_idx] = sortrows(cell2mat(TrialInfo(:,[5 7 9])));

TrialInfo = TrialInfo(sort_idx,:);
TimeVec1024 = Preproc_Data_1.time;

Data_EEG = [Preproc_Data_2.trial; Preproc_Data_1.trial];

fprintf('\n')
fprintf('Artefact Rejection for Subject %s!', Name)
fprintf('\n')
load ArtifactRej
for i = 1:size(trial2rm.(Name),1)
    if(~strcmp(trial2rm.(Name){i,2},'all'))
        Data_EEG(trial2rm.(Name){i,2},trial2rm.(Name){i,1},:) = NaN;
    end
end
curROI(ismember(curROI, cell2mat(trial2rm.(Name)(strcmp(trial2rm.(Name)(:,2),'all'),1)))) = [];

% Fill in NaN Trials
for i = 1:length(Nan_trial_idx)
    Data_EEG = [Data_EEG(1:(Nan_trial_idx(i)-1),:,:); nan(1,size(Data_EEG,2),size(Data_EEG,3)); Data_EEG(Nan_trial_idx(i):end,:,:)];
end

% Sort Data
Data_EEG = Data_EEG(sort_idx,:,:);


%% Test if align
trl_comp = nan(size(TrialInfo,1),1);
for i = 1:size(TrialInfo,1)
    if(TrialInfo{i,1} > 0)
        if(strcmp(TrialInfo{i,6},'Drawing'))
            trl_comp(i,1) = nansum(nansum(squeeze(Preproc_Data_2.trial(cell2mat(Preproc_Data_2.trialinfo(:,2)) == TrialInfo{i,2},:,:) - Data_EEG(i,:,:)),1),2);
        else
            trl_comp(i,1) = nansum(nansum(squeeze(Preproc_Data_1.trial(cell2mat(Preproc_Data_1.trialinfo(:,2)) == TrialInfo{i,2},:,:) - Data_EEG(i,:,:)),1),2);
        end
    end
end

%% Plot features
samplingrate = 1024; %Hz
slide_window_idx = round(samplingrate*(0.020/2));
slide_step_idx = round(samplingrate*0.010);
TimeVec_idx = 1:slide_step_idx:size(Data_EEG,3);
TimeVec = TimeVec1024(TimeVec_idx);

cur_trial1 = zeros(size(Data_EEG,2),1);
cur_trial2 = zeros(size(Data_EEG,2),1);
tp = 40;
i = 1; j = 2;
time_window = max(TimeVec_idx(tp) - slide_window_idx, 1):min(TimeVec_idx(tp) + slide_window_idx, size(Data_EEG,3));
cur_trial1 = squeeze(Data_EEG(i,curROI,time_window));
cur_trial2 = squeeze(Data_EEG(j,curROI,time_window));

figure
for sbp = 1:15
    subplot(5,3,sbp)
    plot(TimeVec1024(time_window), cur_trial1(sbp,:))
    hold on
    plot(TimeVec1024(time_window), cur_trial2(sbp,:))
end

% Euclidian Distance Matrix (Walther et al., 2016)
feature_mat = squeeze(Data_EEG(:,curROI,381));
feature_mat(isnan(feature_mat)) = 0;

c = zeros(1,128);
c(1,[1 5]) = [1 -1]; 

eucl_dist = c*feature_mat*feature_mat'*c';




%% Get Image Data

%TrialInfo = TrialInfo_subs{1,1};
Images = {};
currentdir = pwd;
for i = 1:size(TrialInfo,1)
    
    [~,message,~] = fileattrib(['object_pics\',lower(TrialInfo{i,6}),'_',TrialInfo{i,10},'\*']);
    filenames = strrep({message([message.directory] == 0).Name}',[currentdir,'\'],'');
    filenames = regexprep(filenames,'\','/');
    Images(i,1) = filenames(TrialInfo{i,11});
    Images{i,2} = imread(Images{i,1});
    
end

for i = 1:size(Images,1)
    Images{i,3} = im2double(imresize(Images{i,2},0.03));
end

Test_Mat = ~eye(16).*rand(16);
%Test_Mat = ~eye(128).*rand(128);

col1 = [0.8 0.0 0.0];
col2 = [0.0 0.0 0.8];

imExMDS = Images(1:8:end,3);
imExMDS([1:4 9:12],2) = {uint8(cat(3, 256*col1(1)*ones(70,70), 256*col1(2)*ones(70,70), 256*col1(3)*ones(70,70)))};
imExMDS([5:8 13:16],2) = {uint8(cat(3, 256*col2(1)*ones(70,70), 256*col2(2)*ones(70,70), 256*col2(3)*ones(70,70)))};
plotpictures = false;


for sub = 1:length(Subj_names)
    RSA_Mat_av(sub,:,:,:) = RSA_Mat_subs{sub,2};
    MDS_Mat_av(sub,:,:,:) = MDS_Mat{sub,2};
end

RSA_Mat = squeeze(nanmean(RSA_Mat_av,1));
MDS_Mat = squeeze(nanmean(MDS_Mat_av,1));
TimeVec = RSA_Time.TimeVec;

RSA_Mat = RSA_Mat_av;
MDS_Mat = MDS_Mat_av;

%X = [ones(32,1);2*ones(32,1);ones(32,1);2*ones(32,1)]; %[ones(64,1); 2*ones(64,1)]; %[ones(32,1);2*ones(32,1);ones(32,1);2*ones(32,1)];
X = [ones(4,1);2*ones(4,1);ones(4,1);2*ones(4,1)]; %[ones(8,1); 2*ones(8,1)]; %[ones(4,1);2*ones(4,1);ones(4,1);2*ones(4,1)];
X_names = {'Animate','Inanimate'}; %{'Drawing','Picture'}; %{'Animate','Inanimate'};
%ticks = (8:8:size(RSA_Mat,2)-1)+0.5;
ticks = (4:4:15)+0.5;
chi2_95perc = chi2inv(0.60,2);
t = linspace(0, 2*pi, 50);
mds_dim1 = 1;
mds_dim2 = 2;

h = figure('Pos', [202 69 1552 918]);

% Plot RDM
ax1 = subplot(1,2,1);
pl_im = imagesc(squeeze(RSA_Mat(:,:,1)));
%pl_im = imagesc(Test_Mat);
pbaspect(ax1,[1 1 1])
set(ax1,'XTickLabel',[],'YTickLabel',[],'xtick',ticks,'ytick',ticks,'TickLength',[0.02 0.02],'LineWidth',1.5)
hold on
x_dist = 0.021; count = 0;
for i = 1:8:size(Images,1)-1
    axes('pos',[.122+x_dist*count 0.802 .0355 .0355])
    imshow(Images{i,3})
    axes('pos',[.122+x_dist*count 0.802+0.03 .0355 .0355])
    imshow(Images{i+1,3})
    count = count + 1;
end
x_dist = 0.0354; count = 0;
for i = 1:8:size(Images,1)-1
    axes('pos',[.08 0.76463-x_dist*count .0355 .0355])
    imshow(Images{i+1,3})
    axes('pos',[.08+0.0205 0.76463-x_dist*count .0355 .0355])
    imshow(Images{i,3})
    count = count + 1;
end
text(120,-1030,'Animate','FontSize',14); text(340,-1030,'Inanimate','FontSize',14);
text(590,-1030,'Animate','FontSize',14); text(810,-1030,'Inanimate','FontSize',14);
text(250,-1090,'Drawing','FontSize',14); text(720,-1090,'Picture','FontSize',14);

text(-100,-700,'Animate','FontSize',14,'rotation',90);
text(-100,-460,'Inanimate','FontSize',14,'rotation',90);
text(-100,-240,'Animate','FontSize',14,'rotation',90); 
text(-100,10,'Inanimate','FontSize',14,'rotation',90);
text(-170,-590,'Drawing','FontSize',14,'rotation',90); 
text(-170,-140,'Picture','FontSize',14,'rotation',90);
time_h = text(-170,-1000,sprintf('%3.0f ms',TimeVec(1)*1000),'FontSize',16);

% Plot MDS
ax2 = subplot(1,2,2);
if(plotpictures)
    pl_mds = scatter(MDS_Mat(:,mds_dim1,1),MDS_Mat(:,mds_dim2,1));
    %set(ax2,'Ylim',[-30 30],'Xlim',[-30 30])
    pbaspect(ax2,[1 1 1])
    set(pl_mds,'Marker','none')
    for m = 1:size(MDS_Mat,1)
        imExMDS{m,3} = image('CData', imExMDS{m,2}, 'XData', MDS_Mat(m,mds_dim1,1)+[-size(imExMDS{m,2},2)/50 size(imExMDS{m,2},2)/50], 'YData', MDS_Mat(m,mds_dim2,1)+[size(imExMDS{m,2},1)/50 -size(imExMDS{m,2},1)/50]);
        imExMDS{m,4} = image('CData', imExMDS{m,1}, 'XData', MDS_Mat(m,mds_dim1,1)+[-size(imExMDS{m,1},2)/50 size(imExMDS{m,1},2)/50], 'YData', MDS_Mat(m,mds_dim2,1)+[size(imExMDS{m,1},1)/50 -size(imExMDS{m,1},1)/50]);
    end
else
    pl_mds = gscatter(MDS_Mat(:,mds_dim1,1),MDS_Mat(:,mds_dim2,1),X);
    pbaspect(ax2,[1 1 1])
    hold on
    cent1 = plot(mean(MDS_Mat(X == 1,mds_dim1,1),1),mean(MDS_Mat(X == 1,mds_dim2,1),1),'x','Color',[1 0 0],'MarkerSize',15,'linewidth',3); 
    [V,D]= eig(chi2_95perc*cov(MDS_Mat(X == 1,[mds_dim1 mds_dim2],1)));
    std_err_elipse = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    elipse1 = patch(std_err_elipse(1,:) + mean(MDS_Mat(X == 1,mds_dim1,1),1), std_err_elipse(2,:) + mean(MDS_Mat(X == 1,mds_dim2,1),1), [1 0 0]);
    set(elipse1, 'FaceAlpha',0.2, 'EdgeAlpha',0)
    cent2 = plot(mean(MDS_Mat(X == 2,mds_dim1,1),1),mean(MDS_Mat(X == 2,mds_dim2,1),1),'x','Color',[0 1 1],'MarkerSize',15,'linewidth',3);
    [V,D]= eig(chi2_95perc*cov(MDS_Mat(X == 2,[mds_dim1 mds_dim2],1)));
    std_err_elipse = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    elipse2 = patch(std_err_elipse(1,:) + mean(MDS_Mat(X == 2,mds_dim1,1),1), std_err_elipse(2,:) + mean(MDS_Mat(X == 2,mds_dim2,1),1), [0 1 1]);
    set(elipse2, 'FaceAlpha',0.2, 'EdgeAlpha',0)
    legend(X_names); set(ax2,'Xlim',[-.8 .8],'Ylim',[-.8 .8])
end


scale = 30;
% Loop over Time Points
for i = 2:length(TimeVec)
    set(time_h,'String',sprintf('%3.0f ms',TimeVec(i)*1000))
    set(pl_im,'CData',squeeze(RSA_Mat(:,:,i)))
    if(plotpictures)
        for m = 1:size(MDS_Mat_av,1)
            set(imExMDS{m,3},'XData',MDS_Mat(m,mds_dim1,i)+[-1 1]*(size(imExMDS{m,2},2)/scale), 'YData',MDS_Mat(m,mds_dim2,i)+[1 -1])
            set(imExMDS{m,4},'XData',MDS_Mat(m,mds_dim1,i)+[-1 1]*(size(imExMDS{m,1},2)/scale), 'YData',MDS_Mat(m,mds_dim2,i)+[1 -1])
        end
        yfac = diff(get(ax2,'ylim'))/diff(get(ax2,'xlim'));
        for m = 1:size(MDS_Mat,1)
            set(imExMDS{m,3},'YData',MDS_Mat(m,mds_dim2,i)+[1 -1]*yfac*size(imExMDS{m,2},1)/scale)
            set(imExMDS{m,4},'YData',MDS_Mat(m,mds_dim2,i)+[1 -1]*yfac*size(imExMDS{m,1},1)/scale)
        end
    else
        set(pl_mds(1), 'XData', MDS_Mat(X == 1,mds_dim1,i), 'YData', MDS_Mat(X == 1,mds_dim2,i))
        set(pl_mds(2), 'XData', MDS_Mat(X == 2,mds_dim1,i), 'YData', MDS_Mat(X == 2,mds_dim2,i))
        set(cent1, 'XData', mean(MDS_Mat(X == 1,mds_dim1,i),1), 'YData', mean(MDS_Mat(X == 1,mds_dim2,i),1))
        [V,D]= eig(chi2_95perc*cov(MDS_Mat(X == 1,[mds_dim1 mds_dim2],i)));
        std_err_elipse = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
        set(elipse1, 'XData', std_err_elipse(1,:) + mean(MDS_Mat(X == 1,mds_dim1,i),1), 'YData', std_err_elipse(2,:) + mean(MDS_Mat(X == 1,mds_dim2,i),1))
        set(cent2, 'XData', mean(MDS_Mat(X == 2,mds_dim1,i),1), 'YData', mean(MDS_Mat(X == 2,mds_dim2,i),1))
        [V,D]= eig(chi2_95perc*cov(MDS_Mat(X == 2,[mds_dim1 mds_dim2],i)));
        std_err_elipse = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
        set(elipse2, 'XData', std_err_elipse(1,:) + mean(MDS_Mat(X == 2,mds_dim1,i),1), 'YData',std_err_elipse(2,:) + mean(MDS_Mat(X == 2,mds_dim2,i),1))
    end
    pause(0.5)
end







%% GA Data

% Subject Names
load('RSA_Data_Enc','Subj_names')

% Measure
msr = 1;

% Create Data Struct for Encoding
tmp_strct_enc = load('RSA_Data_Enc');
measures = tmp_strct_enc.RSA_Data_CM.OCC.meas16;
tmp_strct_ret = load('RSA_Data_Ret');
RSA_Data_Enc = [];
RSA_Data_Ret = [];
for sub = 1:length(Subj_names)
    if(sub == 1)
        RSA_Data_Enc.Names   = Subj_names;
        for fn = fieldnames(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC)'
            RSA_Data_Enc.(fn{1}) = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.(fn{1});
        end
        RSA_Data_Enc = rmfield(RSA_Data_Enc, {'Name','TimeVec1024','RSA_full','MDS_full','RSA_16','Encoding_Data','TrialInfo'...
                                      'rsa_dim','curROI','curROI_name'});
        RSA_Data_Enc.OCC_ROI = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.curROI; 
        RSA_Data_Enc.TMP_ROI = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).TMP.curROI;
        RSA_Data_Enc.FRT_ROI = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).FRT.curROI;
        RSA_Data_Enc.PRT_ROI = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).PRT.curROI;

        RSA_Data_Ret.Names   = Subj_names;
        for fn = fieldnames(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC)'
            RSA_Data_Ret.(fn{1}) = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.(fn{1});
        end
        RSA_Data_Ret = rmfield(RSA_Data_Ret, {'Name','TimeVec1024','RSA_full','MDS_full','RSA_16','Retrieval_Data','TrialInfo'...
                                      'rsa_dim','curROI','curROI_name'});
        RSA_Data_Ret.OCC_ROI = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.curROI; 
        RSA_Data_Ret.TMP_ROI = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).TMP.curROI;
        RSA_Data_Ret.FRT_ROI = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).FRT.curROI;
        RSA_Data_Ret.PRT_ROI = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).PRT.curROI;

    end

    if(~isempty(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full))
        RSA_Data_Enc.OCC.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Enc.TMP.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Enc.FRT.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Enc.PRT.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_full{1,msr},[3 1 2]);
    end
    RSA_Data_Enc.OCC.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Enc.TMP.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Enc.FRT.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Enc.PRT.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_16{1,msr},[3 1 2]);

    if(~isempty(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full))
        RSA_Data_Ret.OCC.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Ret.TMP.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Ret.FRT.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Ret.PRT.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_full{1,msr},[3 1 2]);
    end
    RSA_Data_Ret.OCC.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Ret.TMP.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Ret.FRT.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Ret.PRT.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_16{1,msr},[3 1 2]);

    RSA_Data_Enc.Encoding_Data{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.Encoding_Data;
    RSA_Data_Enc.TrialInfo{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;

    RSA_Data_Ret.Retrieval_Data{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.Retrieval_Data;
    RSA_Data_Ret.TrialInfo{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;

end
clear tmp_strct_enc tmp_strct_ret


% ROI Electrode positions
% Biosemi 128 Electrodes System Radial ABC
Electrodes_ROIs

% Hypotheses Matrix
trl_mat = [kron([1;2],ones(64,1)) kron([1;2;1;2],ones(32,1))];
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


RSA_Data = RSA_Data_Enc;

% Get Image Data
TrialInfo = RSA_Data.TrialInfo{1,1};
Images = {};
currentdir = pwd;
for i = 1:size(TrialInfo,1)
    
    [~,message,~] = fileattrib(['object_pics\',lower(TrialInfo{i,6}),'_',TrialInfo{i,10},'\*']);
    filenames = strrep({message([message.directory] == 0).Name}',[currentdir,'\'],'');
    filenames = regexprep(filenames,'\','/');
    Images(i,1) = filenames(TrialInfo{i,11});
    Images{i,2} = imread(Images{i,1});
    
end

for i = 1:size(Images,1)
    Images{i,3} = im2double(imresize(Images{i,2},0.03));
end


ROI = {'OCC','TMP','FRT','PRT'};
Dim_names = {};

[a,b] = unique(TrialInfo(1:64,6),'stable');
Dim_names{1,1} = TrialInfo([b;b+64],6);
Dim_names{1,2} = [8.5; 8.5];
[a,b] = unique(TrialInfo(1:64,8),'stable');
Dim_names{2,1} = TrialInfo([b;b+64],8);
Dim_names{2,2} = [4.5:4:12.5; 4.5:4:12.5];

col1 = [0.8 0.0 0.0];
col2 = [0.0 0.0 0.8];

imExMDS = Images(1:8:end,3);
imExMDS([1:4 9:12],2) = {uint8(cat(3, 256*col1(1)*ones(70,70), 256*col1(2)*ones(70,70), 256*col1(3)*ones(70,70)))};
imExMDS([5:8 13:16],2) = {uint8(cat(3, 256*col2(1)*ones(70,70), 256*col2(2)*ones(70,70), 256*col2(3)*ones(70,70)))};
plotpictures = true;

add_centroid = true;

r = 2;

TimeVec = RSA_Data.TimeVec;
CurRSA = RSA_Data.(ROI{r}).red16_Data;

RSA_Mat = squeeze(mean(CurRSA,1));
tmp = sort(RSA_Mat(:));
rsa_c_limits = [tmp(find(tmp > 0,1,'first'))*0.95 tmp(end)*0.9];

dims = 2;
MDS_Mat = zeros(size(RSA_Mat,2),dims,length(TimeVec));
for tp = 1:length(TimeVec)
    MDS_Mat(:,:,tp) = cmdscale(squeeze(RSA_Mat(tp,:,:)) + squeeze(RSA_Mat(tp,:,:))',dims);
end
mds_x_limits = [min(min(squeeze(MDS_Mat(:,1,:)))) max(max(squeeze(MDS_Mat(:,1,:))))];
mds_y_limits = [min(min(squeeze(MDS_Mat(:,2,:)))) max(max(squeeze(MDS_Mat(:,2,:))))];

chi2_95perc = chi2inv(0.80,2);
t = linspace(0, 2*pi, 50);
mds_dim1 = 1;
mds_dim2 = 2;

Hyp_Mat = cat(3, Semantic_Mat_red16 > 0, Semantic_Mat_red16 < 0);

ts_data = zeros(size(CurRSA,1), size(CurRSA,2), size(Hyp_Mat,3));
for ts = 1:size(Hyp_Mat,3)
    tmp_hyp = Hyp_Mat(:,:,ts);
    ts_data(:,:,ts) = mean(CurRSA(:,:,tmp_hyp(:)),3);
end

X = kron([1;2;1;2], ones(4,1));

posfig = [200 70 1300 900];
h = figure('Pos', posfig);

% Plot RDM
pos1 = [100 340 450 450];
ax1 = axes('Units','pixels','Pos',pos1);
cur_data = squeeze(RSA_Mat(1,:,:));
pl_im = imagesc(squeeze(RSA_Mat(1,:,:)));
clb = colorbar('Units','pixels','Position',[pos1(1)-30 pos1(2) 20 pos1(4)]);
grid on
set(get(gca,'Yruler'),'Minortick',Dim_names{2,2}(1,:))
set(get(gca,'Xruler'),'Minortick',Dim_names{2,2}(1,:))
set(gca,'clim',rsa_c_limits,'xtick',Dim_names{1,2}(1,:),'xticklabels',[],...
        'ytick',Dim_names{1,2}(2,:),'yticklabels',[], 'TickLength',[0 0],'XMinorgrid','on','YMinorgrid','on',...
        'gridcolor','w','gridalpha',.9,'minorgridlinestyle','--','minorgridalpha',.5,'MinorGridColor','w')
count = 0; image_h = [];
for i = 1:8:size(Images,1)-1
    image_h = [image_h, axes('Units','pixels','Pos',[pos1(1)+(pos1(3)/16)*count pos1(2)+pos1(4)+2 pos1(3)/16 pos1(4)/16])];
    imshow(Images{i,3})
    image_h = [image_h, axes('Units','pixels','Pos',[pos1(1)+(pos1(3)/16)*count pos1(2)+pos1(4)+2+pos1(4)/16 pos1(3)/16 pos1(4)/16])];
    imshow(Images{i+1,3})
    image_h = [image_h, axes('Units','pixels','Pos',[pos1(1)+pos1(3)+2  pos1(2)+pos1(4)-pos1(4)/16*(count+1)   pos1(3)/16   pos1(4)/16])];
    imshow(Images{i,3})
    image_h = [image_h, axes('Units','pixels','Pos',[pos1(1)+pos1(3)+2+(pos1(3)/16)  pos1(2)+pos1(4)-pos1(4)/16*(count+1)   pos1(3)/16   pos1(4)/16])];
    imshow(Images{i+1,3})
    count = count + 1;
end
% uicontrol('Style','text','units','pixels','position',[pos1(1)+ 0*(pos1(3)/16) pos1(2)+pos1(4)+8+2*(pos1(4)/16) 4*pos1(3)/16 20],'String','Animate','FontSize',14)
% uicontrol('Style','text','units','pixels','position',[pos1(1)+ 4*(pos1(3)/16) pos1(2)+pos1(4)+8+2*(pos1(4)/16) 4*pos1(3)/16 20],'String','Inanimate','FontSize',14)
% uicontrol('Style','text','units','pixels','position',[pos1(1)+ 8*(pos1(3)/16) pos1(2)+pos1(4)+8+2*(pos1(4)/16) 4*pos1(3)/16 20],'String','Animate','FontSize',14)
% uicontrol('Style','text','units','pixels','position',[pos1(1)+12*(pos1(3)/16) pos1(2)+pos1(4)+8+2*(pos1(4)/16) 4*pos1(3)/16 20],'String','Inanimate','FontSize',14)
% uicontrol('Style','text','units','pixels','position',[pos1(1)+ 0*(pos1(3)/16) pos1(2)+pos1(4)+34+2*(pos1(4)/16) 8*pos1(3)/16 20],'String','Drawing','FontSize',14)
% uicontrol('Style','text','units','pixels','position',[pos1(1)+ 8*(pos1(3)/16) pos1(2)+pos1(4)+34+2*(pos1(4)/16) 8*pos1(3)/16 20],'String','Picture','FontSize',14)
text('parent',ax1,  'String','Animate',   'position',[ 2.5 -2.1 0], 'FontSize',14,  'HorizontalAlignment','center')
text('parent',ax1,  'String','Inanimate', 'position',[ 6.5 -2.1 0], 'FontSize',14,  'HorizontalAlignment','center')
text('parent',ax1,  'String','Animate',   'position',[10.5 -2.1 0], 'FontSize',14,  'HorizontalAlignment','center')
text('parent',ax1,  'String','Inanimate', 'position',[14.5 -2.1 0], 'FontSize',14,  'HorizontalAlignment','center')
text('parent',ax1,  'String','Drawing',   'position',[ 4.5 -3   0], 'FontSize',14,  'HorizontalAlignment','center')
text('parent',ax1,  'String','Picture',   'position',[12.5 -3   0], 'FontSize',14,  'HorizontalAlignment','center')
text('parent',ax1,  'String','Animate',   'position',[19.1  2.5 0], 'FontSize',14,  'HorizontalAlignment','center', 'Rotation',-90)
text('parent',ax1,  'String','Inanimate', 'position',[19.1  6.5 0], 'FontSize',14,  'HorizontalAlignment','center', 'Rotation',-90)
text('parent',ax1,  'String','Animate',   'position',[19.1 10.5 0], 'FontSize',14,  'HorizontalAlignment','center', 'Rotation',-90)
text('parent',ax1,  'String','Inanimate', 'position',[19.1 14.5 0], 'FontSize',14,  'HorizontalAlignment','center', 'Rotation',-90)
text('parent',ax1,  'String','Drawing',   'position',[20    4.5 0], 'FontSize',14,  'HorizontalAlignment','center', 'Rotation',-90)
text('parent',ax1,  'String','Picture',   'position',[20   12.5 0], 'FontSize',14,  'HorizontalAlignment','center', 'Rotation',-90)
time_h = text('parent',ax1,'String',sprintf('%3.0f ms',TimeVec(1)*1000),'position',[8.5 17.1 0],'FontSize',16,'HorizontalAlignment','center');

% Plot MDS
pos2 = [posfig(3)-pos1(3)-pos1(1) pos1(2) pos1(3) pos1(4)];
ax2 = axes('Units','pixels','Pos',pos2);
if(plotpictures)
    scale_mds = 100;
    pl_mds = scatter(MDS_Mat(:,mds_dim1,1).*scale_mds,MDS_Mat(:,mds_dim2,1).*scale_mds);
    %set(ax2,'Ylim',[-30 30],'Xlim',[-30 30])
    pbaspect(ax2,[1 1 1])
    set(pl_mds,'Marker','none')
    for m = 1:size(MDS_Mat,1)
        imExMDS{m,3} = image('CData', imExMDS{m,2}, 'XData', scale_mds.*MDS_Mat(m,mds_dim1,1)+[-size(imExMDS{m,2},2)/50 size(imExMDS{m,2},2)/50], 'YData', scale_mds.*MDS_Mat(m,mds_dim2,1)+[size(imExMDS{m,2},1)/50 -size(imExMDS{m,2},1)/50]);
        imExMDS{m,4} = image('CData', imExMDS{m,1}, 'XData', scale_mds.*MDS_Mat(m,mds_dim1,1)+[-size(imExMDS{m,1},2)/50 size(imExMDS{m,1},2)/50], 'YData', scale_mds.*MDS_Mat(m,mds_dim2,1)+[size(imExMDS{m,1},1)/50 -size(imExMDS{m,1},1)/50]);
    end
else
    scale_mds = 1;
    pl_mds = gscatter(MDS_Mat(:,mds_dim1,1),MDS_Mat(:,mds_dim2,1),X);
    pbaspect(ax2,[1 1 1])
    %legend(X_names); 
end
if(add_centroid)
    hold on
    cent1 = plot(mean(scale_mds.*MDS_Mat(X == 1,mds_dim1,1),1),mean(scale_mds.*MDS_Mat(X == 1,mds_dim2,1),1),'x','Color',[1 0 0],'MarkerSize',15,'linewidth',3); 
    [V,D]= eig(chi2_95perc*cov(scale_mds.*MDS_Mat(X == 1,[mds_dim1 mds_dim2],1)));
    std_err_elipse = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    elipse1 = patch(std_err_elipse(1,:) + mean(scale_mds.*MDS_Mat(X == 1,mds_dim1,1),1), std_err_elipse(2,:) + mean(scale_mds.*MDS_Mat(X == 1,mds_dim2,1),1), [1 0 0]);
    set(elipse1, 'FaceAlpha',0.1, 'EdgeAlpha',0)
    cent2 = plot(mean(scale_mds.*MDS_Mat(X == 2,mds_dim1,1),1),mean(scale_mds.*MDS_Mat(X == 2,mds_dim2,1),1),'x','Color',[0 1 1],'MarkerSize',15,'linewidth',3);
    [V,D]= eig(chi2_95perc*cov(scale_mds.*MDS_Mat(X == 2,[mds_dim1 mds_dim2],1)));
    std_err_elipse = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    elipse2 = patch(std_err_elipse(1,:) + mean(scale_mds.*MDS_Mat(X == 2,mds_dim1,1),1), std_err_elipse(2,:) + mean(scale_mds.*MDS_Mat(X == 2,mds_dim2,1),1), [0 1 1]);
    set(elipse2, 'FaceAlpha',0.1, 'EdgeAlpha',0)
    hold off
end
set(ax2,'Xlim',mds_x_limits.*scale_mds,'Ylim',mds_y_limits.*scale_mds)

% Plot Time Series
pos3 = [pos1(1)+pos1(3)/2 50 pos2(1)+pos2(3)/2 - (pos1(1)+pos1(3)/2) pos1(2)-130];
ax3 = axes('Units','pixels','Pos',pos3);
hold on
dat1 = nanmean(ts_data(:,:,1),1);
dat2 = nanmean(ts_data(:,:,2),1);
SEM1 = nanstd(ts_data(:,:,1),0,1)./sqrt(size(ts_data,1));
SEM2 = nanstd(ts_data(:,:,2),0,1)./sqrt(size(ts_data,1));
fill([TimeVec fliplr(TimeVec)],[dat1 fliplr(dat1 + SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[dat1 fliplr(dat1 - SEM1)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[dat2 fliplr(dat2 + SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
fill([TimeVec fliplr(TimeVec)],[dat2 fliplr(dat2 - SEM2)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
h1 = plot(TimeVec, dat1,'b','linewidth',2);
h2 = plot(TimeVec, dat2,'r','linewidth',2);
xlim([TimeVec(1) TimeVec(end)])
ylimits = get(ax3,'ylim');
time_line = plot([TimeVec(1) TimeVec(1)],[ylimits(1) ylimits(2)],'color',[1 0 1 0.3],'linewidth',2.5);
hold off
set(ax3,'ylim',ylimits)
%lg = legend([h1 h2], {Dat_names{dt(1)},Dat_names{dt(2)}}); legend boxoff; set(lg,'FontSize',14)

set([ax1 ax2 ax3],'units','norm')
set(image_h,'units','norm')
set(clb,'units','norm')

scale = 30;
% Loop over Time Points
for i = 2:length(TimeVec)
    set(time_h,'String',sprintf('%3.0f ms',TimeVec(i)*1000))
    set(time_line,'XData',[TimeVec(i) TimeVec(i)])
    set(pl_im,'CData',squeeze(RSA_Mat(i,:,:)))
    if(plotpictures)
        for m = 1:size(MDS_Mat,1)
            set(imExMDS{m,3},'XData',scale_mds.*MDS_Mat(m,mds_dim1,i)+[-1 1]*(size(imExMDS{m,2},2)/scale), 'YData',scale_mds.*MDS_Mat(m,mds_dim2,i)+[1 -1])
            set(imExMDS{m,4},'XData',scale_mds.*MDS_Mat(m,mds_dim1,i)+[-1 1]*(size(imExMDS{m,1},2)/scale), 'YData',scale_mds.*MDS_Mat(m,mds_dim2,i)+[1 -1])
        end
        yfac = diff(get(ax2,'ylim'))/diff(get(ax2,'xlim'));
        for m = 1:size(MDS_Mat,1)
            set(imExMDS{m,3},'YData',scale_mds.*MDS_Mat(m,mds_dim2,i)+[1 -1]*yfac*size(imExMDS{m,2},1)/scale)
            set(imExMDS{m,4},'YData',scale_mds.*MDS_Mat(m,mds_dim2,i)+[1 -1]*yfac*size(imExMDS{m,1},1)/scale)
        end
    else
        set(pl_mds(1), 'XData', MDS_Mat(X == 1,mds_dim1,i), 'YData', MDS_Mat(X == 1,mds_dim2,i))
        set(pl_mds(2), 'XData', MDS_Mat(X == 2,mds_dim1,i), 'YData', MDS_Mat(X == 2,mds_dim2,i))
    end
    if(add_centroid)
        set(cent1, 'XData', mean(scale_mds.*MDS_Mat(X == 1,mds_dim1,i),1), 'YData', mean(scale_mds.*MDS_Mat(X == 1,mds_dim2,i),1))
        [V,D]= eig(chi2_95perc*cov(scale_mds.*MDS_Mat(X == 1,[mds_dim1 mds_dim2],i)));
        std_err_elipse = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
        set(elipse1, 'XData', std_err_elipse(1,:) + mean(scale_mds.*MDS_Mat(X == 1,mds_dim1,i),1), 'YData', std_err_elipse(2,:) + mean(scale_mds.*MDS_Mat(X == 1,mds_dim2,i),1))
        set(cent2, 'XData', mean(scale_mds.*MDS_Mat(X == 2,mds_dim1,i),1), 'YData', mean(scale_mds.*MDS_Mat(X == 2,mds_dim2,i),1))
        [V,D]= eig(chi2_95perc*cov(scale_mds.*MDS_Mat(X == 2,[mds_dim1 mds_dim2],i)));
        std_err_elipse = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
        set(elipse2, 'XData', std_err_elipse(1,:) + mean(scale_mds.*MDS_Mat(X == 2,mds_dim1,i),1), 'YData',std_err_elipse(2,:) + mean(scale_mds.*MDS_Mat(X == 2,mds_dim2,i),1))
    end
    pause(0.3)
end


