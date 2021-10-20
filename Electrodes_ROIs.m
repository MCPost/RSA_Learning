%% Electrode Names and ROIs

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

% Frontal Cortex ROI
ROI_front = {'D8','D9','D10','D11','D12','D13','D2','D7','D6','D5','D4','D3','C30','C31','C32','C29','C28',...
             'C27','C26','C25','C24','C17','C18','C19','C20','C21','C22','C23','C16','C15','C14','C13','C12',...
             'C11','C2','C8','C9','C10','C4','C3','C7','C6','C5','B27','B28','B29','B30','B31','B32'};
ROI_front_idx = find(cell2mat(cellfun(@(x) any(strcmp(x, ROI_front)), elecs, 'UniformOutput', 0)));

% Parietal Cortex ROI
ROI_pari = {'D24','D25','D26','D27','D28','D17','D16','D32','D31','D30','D29','A7','A6','A12','A11','A10','A9',...
            'A8','A18','A5','A13','A14','A15','A16','A17','A18','A5','A25','A24','A23','A22','A21','A20','A19',...
            'A4','A3','A26','A27','A28','A29','A30','A31','A32','B9','B8','B7','B6','B5','B4','B3','B2','B10',...
            'B10','B11','B12','B13','B14','B15','B16','B17','B18','B19'};
ROI_pari_idx = find(cell2mat(cellfun(@(x) any(strcmp(x, ROI_pari)), elecs, 'UniformOutput', 0)));




cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.center = 'no';
cfg.projection  = 'polar';
cfg.overlap     = 'shift';
cfg.skipscale   = 'yes';
cfg.skipcomnt   = 'yes';
layout_elec = ft_prepare_layout(cfg);
% 
% figure('Pos',[694  464  440  350])
% axes('Units','Pixel','Position',[0  5  435  340]);
% ft_plot_layout(layout_elec,'label','no','box','no','pointsymbol','s','pointcolor','k','pointsize',18)
% h = get(gca,'Child');
% hold on
% plot(h(6).XData(ROI_occ_idx), h(6).YData(ROI_occ_idx),'s','MarkerSize',20,'MarkerFaceColor','k')
% hold off
% set(gca,'Units','norm')
% saveas(gcf,'Results/OCC_ROI.png')
% close(gcf)
% 
% figure('Pos',[694  464  440  350])
% axes('Units','Pixel','Position',[0  5  435  340]);
% ft_plot_layout(layout_elec,'label','no','box','no','pointsymbol','s','pointcolor','k','pointsize',18)
% h = get(gca,'Child');
% hold on
% plot(h(6).XData(ROI_temp_idx), h(6).YData(ROI_temp_idx),'s','MarkerSize',20,'MarkerFaceColor','k')
% hold off
% set(gca,'Units','norm')
% saveas(gcf,'Results/TEMP_ROI.png')
% close(gcf)
% 
% figure('Pos',[694  464  440  350])
% axes('Units','Pixel','Position',[0  5  435  340]);
% ft_plot_layout(layout_elec,'label','no','box','no','pointsymbol','s','pointcolor','k','pointsize',10)
% h = get(gca,'Child');
% hold on
% plot(h(6).XData(ROI_front_idx), h(6).YData(ROI_front_idx),'s','MarkerSize',12,'MarkerFaceColor','k')
% hold off
% set(gca,'Units','norm')
% saveas(gcf,'Results/FRT_ROI.png')
% close(gcf)
% 
% figure('Pos',[694  464  440  350])
% axes('Units','Pixel','Position',[0  5  435  340]);
% ft_plot_layout(layout_elec,'label','no','box','no','pointsymbol','s','pointcolor','k','pointsize',10)
% h = get(gca,'Child');
% hold on
% plot(h(6).XData(ROI_pari_idx), h(6).YData(ROI_pari_idx),'s','MarkerSize',18,'MarkerFaceColor','k')
% hold off
% set(gca,'Units','norm')
% saveas(gcf,'Results/PRT_ROI.png')
% close(gcf)
% 
% 
figure('Pos', [226  57  1185  926])
figpos = get(gcf,'position'); bord_pad = 25;
axes('Units','Pixel','Position',[bord_pad figpos(4)/2+bord_pad figpos(3)/2-2*bord_pad figpos(4)/2-2*bord_pad]);
ft_plot_layout(layout_elec,'label','no','box','no','pointsymbol','s','pointcolor','k','pointsize',10)
h = get(gca,'Child');
hold on
plot(h(6).XData(ROI_occ_idx), h(6).YData(ROI_occ_idx),'s','MarkerSize',12,'MarkerFaceColor','k')
hold off
axes('Units','Pixel','Position',[figpos(3)/2+bord_pad figpos(4)/2+bord_pad figpos(3)/2-2*bord_pad figpos(4)/2-2*bord_pad]);
ft_plot_layout(layout_elec,'label','no','box','no','pointsymbol','s','pointcolor','k','pointsize',10)
h = get(gca,'Child');
hold on
plot(h(6).XData(ROI_temp_idx), h(6).YData(ROI_temp_idx),'s','MarkerSize',12,'MarkerFaceColor','k')
hold off
axes('Units','Pixel','Position',[bord_pad bord_pad figpos(3)/2-2*bord_pad figpos(4)/2-2*bord_pad]);
ft_plot_layout(layout_elec,'label','no','box','no','pointsymbol','s','pointcolor','k','pointsize',10)
h = get(gca,'Child');
hold on
plot(h(6).XData(ROI_front_idx), h(6).YData(ROI_front_idx),'s','MarkerSize',12,'MarkerFaceColor','k')
hold off
axes('Units','Pixel','Position',[figpos(3)/2+bord_pad bord_pad figpos(3)/2-2*bord_pad figpos(4)/2-2*bord_pad]);
ft_plot_layout(layout_elec,'label','no','box','no','pointsymbol','s','pointcolor','k','pointsize',10)
h = get(gca,'Child');
hold on
plot(h(6).XData(ROI_pari_idx), h(6).YData(ROI_pari_idx),'s','MarkerSize',12,'MarkerFaceColor','k')
hold off




