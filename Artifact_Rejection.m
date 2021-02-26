%% Artifact Rejection

trial2rm = [];

%% AG

load Prepro_2_encoding_full_AG.mat r_1_AG r_2_AG

Data_EEG = [r_2_AG.trial; r_1_AG.trial];

% Look at all trials
figure
n = 54;
for i = 1:6
    subplot(3,2,i)
    plot(r_2_AG.time, squeeze(Data_EEG(:,n+i,:)))
    title(num2str(n+i))
end
chan = 125;
idx = find(sum(squeeze(Data_EEG(:,chan,:)) < -50,2) > 0)
figure
plot(r_2_AG.time, squeeze(Data_EEG(idx(1),chan,:)))
if(length(idx) > 1)
    hold on
    for i = 2:length(idx)
        plot(r_2_AG.time, squeeze(Data_EEG(idx(i),chan,:)))
    end
    hold off
end
% Trial (1: 26), (3: 43), (5: 45), (7: 70), (9: 72), (36: 46), (50: 101),
% (55: 10), (57: 110), (58: Total), (60: 64), (63: 94), (81: 18), (84: 63),
% (86: 106), (91: 100), (92: 102), (107: 61), (108: 93), (119: Total),
% (122: 8), (123: 18), (124: 6, 84)

trial2rm.AG =  {1, 3, 5, 7, 9, 36, 50, 55, 57, 58, 60, 63, 81, 84, 86, 91, 92, 107, 108, 119, 122, 123, 124;...
                26, 43, 45, 70, 72, 46, 101, 10, 110, 'all', 64, 94, 18, 63, 106, 100, 102, 61, 93, 'all', 8, 18, [6, 84]}';

%% CEWD

load Prepro_2_encoding_full_CEWD.mat r_1_CEWD r_2_CEWD

Data_EEG = [r_2_CEWD.trial; r_1_CEWD.trial];

% Look at all trials
figure
n = 75;
for i = 1:6
    subplot(3,2,i)
    plot(r_2_CEWD.time, squeeze(Data_EEG(:,n+i,:)))
    title(num2str(n+i))
end
chan = 77;
idx = find(sum(squeeze(Data_EEG(:,chan,:)) < -50,2) > 0)
figure
plot(r_2_CEWD.time, squeeze(Data_EEG(idx(1),chan,:)))
if(length(idx) > 1)
    hold on
    for i = 2:length(idx)
        plot(r_2_CEWD.time, squeeze(Data_EEG(idx(i),chan,:)))
    end
    hold off
end
% Trial (40: 109), (41: 109), (42: 109), (43: 109), (44: 109), (45: 109),
% (46: 109), (47: 68, 109), (48: 109), (59: 83), (74: 11), (77: 115), (96: 56)


trial2rm.CEWD =  {40, 41, 42, 43, 44, 45, 46, 47, 48, 59, 74, 77, 96;...
                  109, 109, 109, 109, 109, 109, 109, 109, 109, 83, 11, 115, 56}';


%% CM

load Prepro_2_encoding_full_CM.mat r_1_CM r_2_CM

Data_EEG = [r_2_CM.trial; r_1_CM.trial];

% Look at all trials
figure('pos',[893 125 1031 791])
n = 126;
for i = 1:6
    subplot(3,2,i)
    plot(r_2_CM.time, squeeze(Data_EEG(:,n+i,:)))
    title(num2str(n+i))
end
chan = 125;
idx = find(sum(squeeze(Data_EEG(:,chan,:)) < -50,2) > 0)
figure
plot(r_2_CM.time, squeeze(Data_EEG(idx(1),chan,:)))
if(length(idx) > 1)
    hold on
    for i = 2:length(idx)
        plot(r_2_CM.time, squeeze(Data_EEG(idx(i),chan,:)))
    end
    hold off
end
% Trial (10: 30,40), (11: 30,40), (13: 43), (14: 43), (15: 43), (16: 43),
% (23: 43), (24: 43), (69: 25), (70: 25), (86: 27,59,84), (87: 47), 
% (89: 19), (128: 91)


trial2rm.CM =  {10, 11, 13, 14, 15, 16, 23, 24, 69, 70, 86, 87, 89, 128;...
                [30,40], [30,40], 43, 43, 43, 43, 43, 43, 25, 25, [27, 59, 84], 47, 19, 91}';


%% DS

load Prepro_2_encoding_full_DS.mat r_1_DS r_2_DS

Data_EEG = [r_2_DS.trial; r_1_DS.trial];

% Look at all trials
figure('pos',[893 125 1031 791])
n = 126;
for i = 1:6
    subplot(3,2,i)
    plot(r_2_DS.time, squeeze(Data_EEG(:,n+i,:)))
    title(num2str(n+i))
end
chan = 125;
idx = find(sum(squeeze(Data_EEG(:,chan,:)) < -50,2) > 0)
figure
plot(r_2_DS.time, squeeze(Data_EEG(idx(1),chan,:)))
if(length(idx) > 1)
    hold on
    for i = 2:length(idx)
        plot(r_2_DS.time, squeeze(Data_EEG(idx(i),chan,:)))
    end
    hold off
end
% Trial (26: 28), (27: 28), (28: 28), (40: 18,28), (41: 28), (56: 86), 
% (57: 86), (58: 86), (73: 86), (120: 16)


trial2rm.DS =  {26, 27, 28, 40, 41, 56, 57, 58, 73, 120;...
                28, 28, 28, [18, 28], 38, 86, 86, 86, 86, 16}';


%% FVM

load Prepro_2_encoding_full_FVM.mat r_1_FVM r_2_FVM

Data_EEG = [r_2_FVM.trial; r_1_FVM.trial];

% Look at all trials
figure('pos',[893 125 1031 791])
n = 126;
for i = 1:6
    subplot(3,2,i)
    plot(r_2_FVM.time, squeeze(Data_EEG(:,n+i,:)))
    title(num2str(n+i))
end
chan = 125;
idx = find(sum(squeeze(Data_EEG(:,chan,:)) < -50,2) > 0)
figure
plot(r_2_FVM.time, squeeze(Data_EEG(idx(1),chan,:)))
if(length(idx) > 1)
    hold on
    for i = 2:length(idx)
        plot(r_2_FVM.time, squeeze(Data_EEG(idx(i),chan,:)))
    end
    hold off
end
% Trial (6: 20), (9: 105), (13: 62), (29: 15), (32: 65), (52: 10,65),
% (81: 27,36), (82: 27,36), (85: 5), (86: 1,2,3,4,5,58,59,60),
% (87: 111), (93: 27,26), (97: 110), (112: 2), (115, 91), (118,38),
% (125: 74), 


trial2rm.FVM =  {6, 9, 13, 29, 32, 52, 81, 82, 85, 86, 87, 93, 97, 112, 115, 118, 125;...
                 20, 105, 62, 15, 65, [10,65], [27,36], [27,36], 5, [1,2,3,4,5,58,59,60], 111, [26,27], 110, 2, 91, 38, 74}';
             


save ArtifactRej trial2rm

            
            