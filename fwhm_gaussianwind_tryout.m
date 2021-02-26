%% FWHM

vec = gausswin(61,2.25);

% Find the half max value.
halfMax = (min(vec) + max(vec)) / 2;
% Find where the vec first drops below half the max.
index1 = find(vec >= halfMax, 1, 'first');
% Find where the vec last rises above half the max.
index2 = find(vec >= halfMax, 1, 'last');
fwhm = index2-index1 + 1; % FWHM in indexes.
1.024*fwhm
plot(vec)


srate = 1024;
time_wind1 = 40;
time_wind2 = 80;
time_scale1 = -time_wind1/2:1:time_wind1/2;
time_scale2 = -time_wind2/2:1:time_wind2/2;

figure('Pos',[16 526 1898 420])
subplot(1,4,1)
data = gausswin(length(time_scale1))'./sum(gausswin(length(time_scale1))); halfMax = (min(data) + max(data)) / 2;
index1 = find(data >= halfMax, 1, 'first');
index2 = find(data >= halfMax, 1, 'last');
plot(time_scale1, data,'k','linewidth',2)
hold on
patch(time_scale1([index1 index2 index2 index1]),[0 0 data([index2 index2])],'m','FaceAlpha',0.1,'EdgeColor','r','LineStyle','--','LineWidth',1);
ylabel('Weight'); xlabel('ms'); title('Time Window 40 ms; alpha = 2.5'); ylim([0 0.06])
text(0,data(index1)*.85,sprintf('FWHM = %2.0f ms',time_scale1(index2) - time_scale1(index1) + 1),'FontSize',13,'HorizontalAlignment','center')

subplot(1,4,2)
data = gausswin(length(time_scale1),1)'./sum(gausswin(length(time_scale1),1)); halfMax = (min(data) + max(data)) / 2;
index1 = find(data >= halfMax, 1, 'first');
index2 = find(data >= halfMax, 1, 'last');
plot(time_scale1, data,'k','linewidth',2)
hold on
patch(time_scale1([index1 index2 index2 index1]),[0 0 data([index2 index2])],'m','FaceAlpha',0.1,'EdgeColor','r','LineStyle','--','LineWidth',1);
xlabel('ms'); title('Time Window 40 ms; alpha = 1'); ylim([0 0.06])
text(0,data(index1)*.85,sprintf('FWHM = %2.0f ms',time_scale1(index2) - time_scale1(index1) + 1),'FontSize',13,'HorizontalAlignment','center')

subplot(1,4,3)
data = gausswin(length(time_scale2))'./sum(gausswin(length(time_scale2))); halfMax = (min(data) + max(data)) / 2;
index1 = find(data >= halfMax, 1, 'first');
index2 = find(data >= halfMax, 1, 'last');
plot(time_scale2, data,'k','linewidth',2)
hold on
patch(time_scale2([index1 index2 index2 index1]),[0 0 data([index2 index2])],'m','FaceAlpha',0.1,'EdgeColor','r','LineStyle','--','LineWidth',1);
xlabel('ms'); title('Time Window 80 ms; alpha = 2.5'); ylim([0 0.03])
text(0,data(index1)*.85,sprintf('FWHM = %2.0f ms',time_scale2(index2) - time_scale2(index1) + 1),'FontSize',13,'HorizontalAlignment','center')

subplot(1,4,4)
data = gausswin(length(time_scale2),1)'./sum(gausswin(length(time_scale2),1)); halfMax = (min(data) + max(data)) / 2;
index1 = find(data >= halfMax, 1, 'first');
index2 = find(data >= halfMax, 1, 'last');
plot(time_scale2, data,'k','linewidth',2)
hold on
patch(time_scale2([index1 index2 index2 index1]),[0 0 data([index2 index2])],'m','FaceAlpha',0.1,'EdgeColor','r','LineStyle','--','LineWidth',1);
xlabel('ms'); title('Time Window 80 ms; alpha = 1'); ylim([0 0.03])
text(0,data(index1)*.85,sprintf('FWHM = %2.0f ms',time_scale2(index2) - time_scale2(index1) + 1),'FontSize',13,'HorizontalAlignment','center')






