%% Compare RSA Matrices Simulation

N = 16;

figure('Pos', [326 71 1291 912])

mu = [-5 -5; 5 5];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,1)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-10 10], 'ylim',[-10 10])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-1 -1; 1 1];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,2)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [0 0; 0 0];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,3)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,2)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-1 1; 1 -1];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,4)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-5 5; 5 -5];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,5)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-10 10], 'ylim',[-10 10])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-5 -5; 5 5];
sigma = cat(3,[1 .9; .9 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,6)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-10 10], 'ylim',[-10 10])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-1 -1; 1 1];
sigma = cat(3,[1 .9; .9 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,7)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [0 0; 0 0];
sigma = cat(3,[1 .9; .9 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,8)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-1 1; 1 -1];
sigma = cat(3,[1 .9; .9 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,9)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-5 5; 5 -5];
sigma = cat(3,[1 .9; .9 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,10)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-10 10], 'ylim',[-10 10])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-5 -5; 5 5];
sigma = cat(3,[1 .9; .9 1], [1 -.9; -.9 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,11)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-10 10], 'ylim',[-10 10])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,2)), max(R(:,2)),N),B(2) + linspace(min(R(:,2)), max(R(:,2)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,2)), max(R(1:N,2)),N),B(2) + linspace(min(R(1:N,2)), max(R(1:N,2)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,2)), max(R(N+1:end,2)),N),B(2) + linspace(min(R(N+1:end,2)), max(R(N+1:end,2)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-1 -1; 1 1];
sigma = cat(3,[1 .9; .9 1], [1 -.9; -.9 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,12)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [0 0; 0 0];
sigma = cat(3,[1 .9; .9 1], [1 -.9; -.9 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,13)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-1 1; 1 -1];
sigma = cat(3,[1 .9; .9 1], [1 -.9; -.9 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,14)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-5 5; 5 -5];
sigma = cat(3,[1 .9; .9 1], [1 -.9; -.9 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,15)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-10 10], 'ylim',[-10 10])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-5 -5; 5 5];
sigma = cat(3,[1 -.9; -.9 1], [1 -.9; -.9 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,16)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-10 10], 'ylim',[-10 10])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,2)), max(R(:,2)),N),B(2) + linspace(min(R(:,2)), max(R(:,2)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,2)), max(R(1:N,2)),N),B(2) + linspace(min(R(1:N,2)), max(R(1:N,2)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,2)), max(R(N+1:end,2)),N),B(2) + linspace(min(R(N+1:end,2)), max(R(N+1:end,2)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-1 -1; 1 1];
sigma = cat(3,[1 -.9; -.9 1], [1 -.9; -.9 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,17)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [0 0; 0 0];
sigma = cat(3,[1 -.9; -.9 1], [1 -.9; -.9 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,18)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-1 1; 1 -1];
sigma = cat(3,[1 -.9; -.9 1], [1 -.9; -.9 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,19)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-5 5], 'ylim',[-5 5])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))

mu = [-5 5; 5 -5];
sigma = cat(3,[1 -.9; -.9 1], [1 -.9; -.9 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

subplot(4,5,20)
plot(R(:,1),R(:,2),'+')
set(gca,'xlim',[-10 10], 'ylim',[-10 10])
hold on
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off
title(sprintf('V_{corr} = %1.4f', (atanh(corr(R(1:N,1), R(1:N,2))) - atanh(corr(R(N+1:end,1), R(N+1:end,2))))*(atanh(corr(R(:,1), R(:,2))))))



%% Matrix 


N = 16;

%figure('Pos', [326 71 1291 912])

mu = [0 0; 0 0];
sigma = cat(3,[1 0.99; 0.99 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),0.5*N*(0.5*N -1))];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),0.5*N*(0.5*N -1))];

fast_corr(tiedrank_(R(:,1),1), tiedrank_(R(:,2),1));
ICC(tiedrank_(R,1), 'A-k');
ICC(tiedrank_(R,1), 'C-k');


corrs1 = -0.9:0.05:0.9;
corrs2 = -0.9:0.05:0.9;

mu = [0 0; 0 0];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
Cor_Mat = zeros(length(corrs1),length(corrs2));
ICC_C_Mat = zeros(length(corrs1),length(corrs2));
ICC_A_Mat = zeros(length(corrs1),length(corrs2));
for c1 = 1:length(corrs1)
    for c2 = 1:length(corrs2)
        sigma(1,2,1) = corrs1(c1); sigma(2,1,1) = corrs1(c1);
        sigma(1,2,2) = corrs2(c2); sigma(2,1,2) = corrs2(c2);
        R = [];
        R = [R; mvnrnd(mu(1,:),sigma(:,:,1),0.5*N*(0.5*N -1))];
        R = [R; mvnrnd(mu(2,:),sigma(:,:,2),0.5*N*(0.5*N -1))];

        Cor_Mat(c1,c2) = fast_corr(tiedrank_(R(:,1),1), tiedrank_(R(:,2),1));
        %ICC_A_Mat(c1,c2) = ICC([R(1:56,1:2) R(57:112,1:2)], 'A-1');
        %ICC_C_Mat(c1,c2) = ICC([R(1:56,1:2) R(57:112,1:2)], 'C-1');
        tr_data = tiedrank_(R,1);
        ICC_A_Mat(c1,c2) = ICC([tr_data(1:56,1:2) tr_data(57:112,1:2)], 'A-k');
        ICC_C_Mat(c1,c2) = ICC([tr_data(1:56,1:2) tr_data(57:112,1:2)], 'C-k');
    end
end

figure
subplot(1,3,1)
surf(corrs1, corrs2, Cor_Mat)
title('Correlation'), xlabel('x'), ylabel('y')
subplot(1,3,2)
surf(corrs1, corrs2, ICC_A_Mat)
%set(gca,'zlim',[-3 1.2])
title('ICC Agreement'), xlabel('x'), ylabel('y')
subplot(1,3,3)
surf(corrs1, corrs2, ICC_C_Mat)
%set(gca,'zlim',[-3 1.2])
title('ICC Consistency'), xlabel('x'), ylabel('y')

figure
surf(corrs1, corrs2, ICC_A_Mat./ICC_C_Mat)
title('ICC Agreement/Consistency'), xlabel('x'), ylabel('y')


mus1 = -3:0.2:3;
mus2 = -3:0.2:3;

mu = [0 0; 0 0];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
Cor_Mat = zeros(length(mus1),length(mus2));
ICC_C_Mat = zeros(length(mus1),length(mus2));
ICC_A_Mat = zeros(length(mus1),length(mus2));
for m1 = 1:length(mus1)
    for m2 = 1:length(mus2) 
        mu(2,:) = [mus1(m1) mus2(m2)];
        R = [];
        R = [R; mvnrnd(mu(1,:),sigma(:,:,1),0.5*N*(0.5*N -1))];
        R = [R; mvnrnd(mu(2,:),sigma(:,:,2),0.5*N*(0.5*N -1))];

        Cor_Mat(m1,m2) = fast_corr(tiedrank_(R(:,1),1), tiedrank_(R(:,2),1));
        %ICC_A_Mat(m1,m2) = ICC([R(1:56,1:2) R(57:112,1:2)], 'A-1');
        %ICC_C_Mat(m1,m2) = ICC([R(1:56,1:2) R(57:112,1:2)], 'C-1');
        tr_data = tiedrank_(R,1);
        ICC_A_Mat(m1,m2) = ICC([tr_data(1:56,1:2) tr_data(57:112,1:2)], 'A-k');
        ICC_C_Mat(m1,m2) = ICC([tr_data(1:56,1:2) tr_data(57:112,1:2)], 'C-k');
    end
end

figure
subplot(1,3,1)
surf(mus1, mus2, Cor_Mat)
title('Correlation'), xlabel('x'), ylabel('y')
subplot(1,3,2)
surf(mus1, mus2, ICC_A_Mat)
%set(gca,'zlim',[-3 1.2])
title('ICC Agreement'), xlabel('x'), ylabel('y')
subplot(1,3,3)
surf(mus1, mus2, ICC_C_Mat)
%set(gca,'zlim',[-3 1.2])
title('ICC Consistency'), xlabel('x'), ylabel('y')


figure
surf(mus1, mus2, ICC_A_Mat./ICC_C_Mat)
title('ICC Agreement/Consistency'), xlabel('x'), ylabel('y')



mus1 = -2:1:2;
mus2 = -2:1:2;

mu = [0 0; 0 0];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
Cor_Cell = cell(length(mus1),length(mus2));
ICC_C_Cell = cell(length(mus1),length(mus2));
ICC_A_Cell = cell(length(mus1),length(mus2));
for m1 = 1:length(mus1)
    for m2 = 1:length(mus2) 

        corrs1 = -0.9:0.1:0.9;
        corrs2 = -0.9:0.1:0.9;

        Cor_Cell{m1,m2} = zeros(length(corrs1),length(corrs2));
        ICC_C_Cell{m1,m2} = zeros(length(corrs1),length(corrs2));
        ICC_A_Cell{m1,m2} = zeros(length(corrs1),length(corrs2));
        for c1 = 1:length(corrs1)
            for c2 = 1:length(corrs2)
                mu(2,:) = [mus1(m1) mus2(m2)];
                sigma(1,2,1) = corrs1(c1); sigma(2,1,1) = corrs1(c1);
                sigma(1,2,2) = corrs2(c2); sigma(2,1,2) = corrs2(c2);
                R = [];
                R = [R; mvnrnd(mu(1,:),sigma(:,:,1),0.5*N*(0.5*N -1))];
                R = [R; mvnrnd(mu(2,:),sigma(:,:,2),0.5*N*(0.5*N -1))];

                Cor_Cell{m1,m2}(c1,c2) = fast_corr(tiedrank_(R(:,1),1), tiedrank_(R(:,2),1));
                ICC_C_Cell{m1,m2}(c1,c2) = ICC([R(1:56,1:2) R(57:112,1:2)], 'A-1');
                ICC_A_Cell{m1,m2}(c1,c2) = ICC([R(1:56,1:2) R(57:112,1:2)], 'C-1');
                %tr_data = tiedrank_(R,1);
                %ICC_C_Cell{m1,m2}(c1,c2) = ICC([tr_data(1:56,1:2) tr_data(57:112,1:2)], 'A-1');
                %ICC_A_Cell{m1,m2}(c1,c2) = ICC([tr_data(1:56,1:2) tr_data(57:112,1:2)], 'C-1');
            end
        end
        
    end
end


figure
Data = ICC_C_Cell';
tmp_max = cell2mat(cellfun(@(x) max(x(:)), Data,'UniformOutput',0));
tmp_min = cell2mat(cellfun(@(x) min(x(:)), Data,'UniformOutput',0));
c_limit = [min(tmp_min(:)) max(tmp_max(:))];
for sbp = 1:25
    subplot(5,5,sbp)
    imagesc(Data{sbp})
    set(gca,'clim',c_limit, 'xtick',[], 'ytick',[])
    axis square
end 

figure
Data1 = ICC_C_Cell';
Data2 = ICC_A_Cell';
tmp_max = cell2mat(cellfun(@(x,y) max(x(:)./y(:)), Data1, Data2,'UniformOutput',0));
tmp_min = cell2mat(cellfun(@(x,y) min(x(:)./y(:)), Data1, Data2,'UniformOutput',0));
c_limit = [min(tmp_min(:)) max(tmp_max(:))];
for sbp = 1:25
    subplot(5,5,sbp)
    cur_data = Data1{sbp}./Data2{sbp};
    imagesc(cur_data)
    set(gca,'clim',c_limit,'xtick',[], 'ytick',[])
    axis square
end 


N = 16;

mus1 = [-3 -1 0 1 3]; %-2:1:2;
mus2 = [3 1 0 -1 -3]; %2:-1:-2;

mu = [0 0; 0 0];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
NMeas_Cell = cell(length(mus2),length(mus1));
for m2 = 1:length(mus2)
    for m1 = 1:length(mus1) 

        corrs1 = -0.9:0.1:0.9;
        corrs2 = 0.9:-0.1:-0.9;

        NMeas_Cell{m2,m1} = zeros(length(corrs2),length(corrs1));
        for c2 = 1:length(corrs2)
            for c1 = 1:length(corrs1)
                mu(2,:) = [mus1(m1) mus2(m2)];
                sigma(1,2,1) = corrs1(c1); sigma(2,1,1) = corrs1(c1);
                sigma(1,2,2) = corrs2(c2); sigma(2,1,2) = corrs2(c2);
                R = [];
                R = [R; mvnrnd(mu(1,:),sigma(:,:,1),0.5*N*(0.5*N -1))];
                R = [R; mvnrnd(mu(2,:),sigma(:,:,2),0.5*N*(0.5*N -1))];
                
                cur_data_all = tiedrank_(R,1);
                cor_all = fast_corr(cur_data_all(:,1), cur_data_all(:,2));
                cur_data = tiedrank_(R(1:56,:),1);
                cor_wi  = fast_corr(cur_data(:,1), cur_data(:,2));
                cur_data = tiedrank_(R(57:112,:),1);
                cor_bt  = fast_corr(cur_data(:,1), cur_data(:,2));
                %NMeas_Cell{m2,m1}(c2,c1) = (atanh(cor_wi) - atanh(cor_bt))*atanh(cor_all);
                %[~,p_val] = ttest2(reshape(R(1:56,:),56*2,1), reshape(R(56+1:end,:),56*2,1));
                %NMeas_Cell{m2,m1}(c2,c1) = (abs(norminv(p_val,0,1))/abs(atanh(cor_wi) - atanh(cor_bt)))*atanh(cor_all);
                
                %X = [ones(N*(0.5*N -1),1), cur_data_all(:,2)];
                %%X = [ones(N*(0.5*N -1),1), kron([1;0], ones(0.5*N*(0.5*N -1),1))];
                %%X = R_tr(:,2);
                %y = cur_data_all(:,1);
                %SS_tot = (y - mean(y))' * (y - mean(y));
                %SS_res1 = y'*(eye(N*(0.5*N -1)) - X*((X'*X)\X'))*y;
                %%Rsq1 = 1 - SS_res1/SS_tot;
                %X = [ones(N*(0.5*N -1),1), cur_data_all(:,2), kron([1;0], ones(0.5*N*(0.5*N -1),1))];
                %%X = [ones(N*(0.5*N -1),1), cur_data_all(:,2), kron([1;0], ones(0.5*N*(0.5*N -1),1)).*cur_data_all(:,2)];
                %%X = [R_tr(:,2), kron([1;0], ones(0.5*N*(0.5*N -1),1)), kron([0;1], ones(0.5*N*(0.5*N -1),1))];
                %SS_res2 = y'*(eye(N*(0.5*N -1)) - X*((X'*X)\X'))*y;
                %%Rsq2 = 1 - SS_res2/SS_tot;
                %F = (SS_res1 - SS_res2)/(SS_res2 / (N*(0.5*N -1) - 3));
                %NMeas_Cell{m2,m1}(c2,c1) = abs(f2z_bloc(F,1,N*(0.5*N -1) - 3))*atanh(cor_all);
                
                %y = (cur_data_all(:,1) - mean(cur_data_all(:))) .* (cur_data_all(:,2) - mean(cur_data_all(:)));
                y = tiedrank_(R(:,1) .* R(:,2),1);
                X = [kron([1;0], ones(0.5*N*(0.5*N -1),1)) kron([0;1], ones(0.5*N*(0.5*N -1),1))];
                SS_tot = (y - mean(y))' * (y - mean(y));
                SS_res = (y - X*((X'*X)\X'*y))' * (y - X*((X'*X)\X'*y));
                NMeas_Cell{m2,m1}(c2,c1) = atanh(sqrt(max(1 - SS_res/SS_tot,0))) * atanh(cor_all);
                
                %[~,~,stats] = manova1(cur_data_all,[repmat({'sub1'},0.5*N*(0.5*N -1),1); repmat({'sub2'},0.5*N*(0.5*N -1),1)]);
                %NMeas_Cell{m2,m1}(c2,c1) = atanh(sqrt(1 - min(stats.lambda,1))) * atanh(cor_all);
                
                %X = [kron([1;0], ones(0.5*N*(0.5*N -1),1)) kron([0;1], ones(0.5*N*(0.5*N -1),1))];
                
                %dat1 = R(logical(X(:,1)),:); dat2 = R(logical(X(:,2)),:);
                %y1_m = mean(dat1,1);
                %y2_m = mean(dat2,1);
                %S1 = cov(dat1);
                %S2 = cov(dat2);
                %S = ((0.5*N*(0.5*N -1)-1)*S1 + (0.5*N*(0.5*N -1)-1)*S2)/(N*(0.5*N -1)-2);
                
                %T2 = ((y1_m - y2_m)/inv(S)*(y1_m - y2_m)')*((0.5*N*(0.5*N -1)*0.5*N*(0.5*N -1))/(N*(0.5*N -1)));
                %F = T2*((N*(0.5*N -1) - 3) / (2*(N*(0.5*N -1) - 2)));
                %NMeas_Cell{m2,m1}(c2,c1) = f2z_bloc(F,1,N*(0.5*N -1) - 3);%*atanh(cor_all);
                
            end
        end
        
    end
end
prctile(NMeas_Cell{3,3}(:),[1 99])

figure('Pos',[164 67 1504 926])
Data = NMeas_Cell';
tmp_max = cell2mat(cellfun(@(x) max(x(:)), Data,'UniformOutput',0));
tmp_min = cell2mat(cellfun(@(x) min(x(:)), Data,'UniformOutput',0));
c_limit = [-0.3 0.3]; %prctile(NMeas_Cell{3,3}(:),[1 99]);%[min(tmp_min(:))+0.9 max(tmp_max(:))-1];
phase_shift = [1.2 0.6 0 -0.6 -1.2];
for sbp = 1:25
    h1 = subplot(5,5,sbp);
    imagesc(Data{sbp});
    set(gca,'clim',c_limit, 'xtick',[], 'ytick',[])
    axis square
    aspect = get(h1,'PlotBoxAspectRatio');
    set(h1,'Units','pixels');
    pos = get(h1,'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(h1,'Position',pos);
    
    h2 = axes('Units','pixels','pos',[pos(1) pos(2)-25 pos(3) 20],'xlim',[-15 155],'ylim',[-15 15],'visible','off');
    hold on
    for elp = 1:5
        t = linspace(0,2*pi) ;
        a = 15; b = 15;
        x = a*cos(t+phase_shift(elp)) ;
        y = b*sin(t) ;
        plot(x+35*(elp-1),y,'k','linewidth',2)
    end
    axis equal
    box off
    h3 = axes('Units','pixels','pos',[pos(1)-25 pos(2) 20 pos(4)],'xlim',[-15 15],'ylim',[-15 155],'visible','off');
    hold on
    for elp = 1:5
        t = linspace(0,2*pi) ;
        a = 15; b = 15;
        x = a*cos(t+phase_shift(elp)) ;
        y = b*sin(t) ;
        plot(x,y+35*(elp-1),'k','linewidth',2)
    end
    axis equal
    box off
    
    set(h1,'Units','normalized')
    set(h2,'Units','normalized')
    set(h3,'Units','normalized')
    
end 


%% Plot Correlations

mus = [kron(ones(5,1),[-3 -1 0 1 3]') kron([-3 -1 0 1 3]',ones(5,1))];

figure('Pos', [326 71 1291 912])

mu = [0 0; 0 0];
sigma = cat(3,[1 .9; .9 1], [1 .9; .9 1]);
for sbp = 1:25
    mu(2,:) = [mus(sbp) mus(sbp)];
    R = [];
    R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
    R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];
    R_tr = tiedrank_(R,1); R_tr_g1 = tiedrank_(R(1:N,:),1); R_tr_g2 = tiedrank_(R(N+1:end,:),1);
    h1 = subplot(5,5,sbp);
    plot(R_tr(1:N,1),R_tr(1:N,2),'+r')
    hold on
    plot(R_tr(N+1:end,1),R_tr(N+1:end,2),'+g')
    set(gca,'xlim',[-5 40], 'ylim',[-5 40], 'xtick', [],'ytick',[]) %set(gca,'xlim',[-5 5], 'ylim',[-5 5], 'xtick', [],'ytick',[])
    B = polyfit(R_tr(:,1), R_tr(:,2),1);
    ls1 = plot(linspace(min(R_tr(:,1)), max(R_tr(:,1)),N),B(2) + linspace(min(R_tr(:,1)), max(R_tr(:,1)),N)*B(1),'--m','linewidth',1.5);
    text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R_tr(:,1), R_tr(:,2))),'0.','.'),'Color','m')
    B = polyfit(R_tr(1:N,1), R_tr(1:N,2),1);
    ls2 = plot(linspace(min(R_tr(1:N,1)), max(R_tr(1:N,1)),N),B(2) + linspace(min(R_tr(1:N,1)), max(R_tr(1:N,1)),N)*B(1),'--r','linewidth',1.5);
    text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R_tr(1:N,1), R_tr(1:N,2))),'0.','.'),'Color','r')
    B = polyfit(R_tr(N+1:end,1), R_tr(N+1:end,2),1);
    ls3 = plot(linspace(min(R_tr(N+1:end,1)), max(R_tr(N+1:end,1)),N),B(2) + linspace(min(R_tr(N+1:end,1)), max(R_tr(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
    text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R_tr(N+1:end,1), R_tr(N+1:end,2))),'0.','.'),'Color','g')
    hold off
    
    title(sprintf('V_{corr} = %1.4f', (atanh(corr(R_tr_g1(:,1), R_tr_g1(:,2))) - atanh(corr(R_tr_g2(:,1), R_tr_g2(:,2))))*(atanh(corr(R_tr(:,1), R_tr(:,2))))),'fontsize',8)
    axis square
end


%% New Idea

N = 56;

mus = [kron(ones(5,1),([-3 -1 0 1 3])') kron(([3 1 0 -1 -3])',ones(5,1))];

figure('Pos', [326 71 1291 912])

mu = [0 0; 0 0];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
for sbp = 1:25
    mu(2,:) = [mus(sbp,1) mus(sbp,2)];
    R = [];
    R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
    R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];
    R_tr = tiedrank_(R,1); R_tr_g1 = tiedrank_(R(1:N,:),1); R_tr_g2 = tiedrank_(R(N+1:end,:),1);
    %R_tr = R;
    h1 = subplot(5,5,sbp);
    plot(R_tr(1:N,1),R_tr(1:N,2),'+r')
    hold on
    plot(R_tr(N+1:end,1),R_tr(N+1:end,2),'+g')
    set(gca,'xlim',[min(R_tr(:))-20 max(R_tr(:))+13], 'ylim',[min(R_tr(:))-20 max(R_tr(:))+13], 'xtick', [],'ytick',[]) %set(gca,'xlim',[-5 5], 'ylim',[-5 5], 'xtick', [],'ytick',[])
    B = polyfit(R_tr(:,1), R_tr(:,2),1);
    ls1 = plot(linspace(min(R_tr(:,1)), max(R_tr(:,1)),N),B(2) + linspace(min(R_tr(:,1)), max(R_tr(:,1)),N)*B(1),'--m','linewidth',1.5);
    text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R_tr(:,1), R_tr(:,2))),'0.','.'),'Color','m')
    B = polyfit(R_tr(1:N,1), R_tr(1:N,2),1);
    ls2 = plot(linspace(min(R_tr(1:N,1)), max(R_tr(1:N,1)),N),B(2) + linspace(min(R_tr(1:N,1)), max(R_tr(1:N,1)),N)*B(1),'--r','linewidth',1.5);
    text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R_tr(1:N,1), R_tr(1:N,2))),'0.','.'),'Color','r')
    B = polyfit(R_tr(N+1:end,1), R_tr(N+1:end,2),1);
    ls3 = plot(linspace(min(R_tr(N+1:end,1)), max(R_tr(N+1:end,1)),N),B(2) + linspace(min(R_tr(N+1:end,1)), max(R_tr(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
    text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R_tr(N+1:end,1), R_tr(N+1:end,2))),'0.','.'),'Color','g')
    hold off
    
    %X = [ones(2*N,1), R_tr(:,2)];
    %%X = [ones(2*N,1), kron([1;0], ones(N,1))];
    %%X = R_tr(:,2);
    %y = R_tr(:,1);
    %SS_tot = (y - mean(y))' * (y - mean(y));
    %SS_res1 = (y - X*((X'*X)\X'*y))' * (y - X*((X'*X)\X'*y));
    %%Rsq1 = 1 - SS_res1/SS_tot;
    %X = [ones(2*N,1), R_tr(:,2), kron([1;0], ones(N,1))];
    %%X = [ones(2*N,1), R_tr(:,2), kron([1;0], ones(N,1)).*R_tr(:,2)];
    %%X = [R_tr(:,2), kron([1;0], ones(16,1)), kron([0;1], ones(16,1))];
    %SS_res2 = (y - X*((X'*X)\X'*y))' * (y - X*((X'*X)\X'*y));
    %%Rsq2 = 1 - SS_res2/SS_tot;
    %F = (SS_res1 - SS_res2)/(SS_res2 / (2*N - 3));
    %p = fcdf(F,1,2*N - 3,'upper');
    %meas = abs(norminv(p))*atanh(corr(R_tr(:,1), R_tr(:,2)));
    
    y = tiedrank_(R(:,1) .* R(:,2),1);
    X = [kron([1;0], ones(N,1)), kron([0;1], ones(N,1))];
    SS_tot = (y - mean(y))' * (y - mean(y));
    SS_res = (y - X*((X'*X)\X'*y))' * (y - X*((X'*X)\X'*y));
    meas = atanh(sqrt(1 - SS_res/SS_tot)) * atanh(corr(R_tr(:,1), R_tr(:,2)));
    
    %[~,p_val] = ttest2(reshape(R(1:N,:),N*2,1), reshape(R(N+1:end,:),N*2,1));
    %meas = (atanh(corr(R_tr_g1(:,1), R_tr_g1(:,2))) - atanh(corr(R_tr_g2(:,1), R_tr_g2(:,2))))*(atanh(corr(R_tr(:,1), R_tr(:,2))))*norminv(1-p_val,0,1);
    
    %X = [kron([1;0], ones(N,1)) kron([0;1], ones(N,1))];

    %dat1 = R(logical(X(:,1)),:); dat2 = R(logical(X(:,2)),:);
    %y1_m = mean(dat1,1);
    %y2_m = mean(dat2,1);
    %S1 = cov(dat1);
    %S2 = cov(dat2);
    %S = ((0.5*N*(0.5*N -1)-1)*S1 + (0.5*N*(0.5*N -1)-1)*S2)/(N*(0.5*N -1)-2);

    %T2 = ((y1_m - y2_m)/inv(S)*(y1_m - y2_m)')*((0.5*N*(0.5*N -1)*0.5*N*(0.5*N -1))/(N*(0.5*N -1)));
    %meas = T2;
    title(sprintf('V_{corr} = %1.4f', meas),'fontsize',8)
    axis square
end



%% New New Crossproduct Idea

N = 56;

mus = [kron(ones(5,1),([-3 -1 0 1 3])') kron(([3 1 0 -1 -3])',ones(5,1))];

figure('Pos', [326 71 1291 912])

mu = [0 0; 0 0];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
for sbp = 1:25
    mu(2,:) = [mus(sbp,1) mus(sbp,2)];
    R = [];
    R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
    R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];
    R_tr = tiedrank_(R,1); 
    dat = tiedrank_(R(:,1).*R(:,2),1); % R_tr(:,1).*R_tr(:,2); R(:,1).*R(:,2);
    h1 = subplot(5,5,sbp);
    hold on
    histogram(dat(1:N),15)
    histogram(dat(N+1:end),15)
    hold off
    
    y = dat;
    X = [kron([1;0], ones(N,1)), kron([0;1], ones(N,1))];
    SS_tot = (y - mean(y))' * (y - mean(y));
    SS_res = (y - X*((X'*X)\X'*y))' * (y - X*((X'*X)\X'*y));
    meas = atanh(sqrt(1 - SS_res/SS_tot)) * atanh(corr(R_tr(:,1), R_tr(:,2)));
    
    title(sprintf('V_{corr} = %1.4f', meas),'fontsize',8)
    set(gca, 'xlim', [-2 130], 'ylim', [0 15])
    axis square
end




%% Try out Real Data

msr = 1;

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

% Load Data 
% Subject Names
load('RSA_Data_Enc','Subj_names')

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
        RSA_Data_Enc.CNT_ROI = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).CNT.curROI;
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
        RSA_Data_Ret.CNT_ROI = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).CNT.curROI;
        RSA_Data_Ret.PRT_ROI = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).PRT.curROI;

    end

    if(~isempty(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full))
        RSA_Data_Enc.OCC.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Enc.TMP.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Enc.FRT.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Enc.CNT.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).CNT.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Enc.PRT.full_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_full{1,msr},[3 1 2]);
    end
    RSA_Data_Enc.OCC.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Enc.TMP.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Enc.FRT.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Enc.CNT.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).CNT.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Enc.PRT.red16_Data(sub,:,:,:) = permute(tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_16{1,msr},[3 1 2]);

    if(~isempty(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full))
        RSA_Data_Ret.OCC.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Ret.TMP.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Ret.FRT.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Ret.CNT.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).CNT.RSA_full{1,msr},[3 1 2]);
        RSA_Data_Ret.PRT.full_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_full{1,msr},[3 1 2]);
    end
    RSA_Data_Ret.OCC.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Ret.TMP.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).TMP.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Ret.FRT.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).FRT.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Ret.CNT.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).CNT.RSA_16{1,msr},[3 1 2]);
    RSA_Data_Ret.PRT.red16_Data(sub,:,:,:) = permute(tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).PRT.RSA_16{1,msr},[3 1 2]);

    RSA_Data_Enc.Encoding_Data{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.Encoding_Data;
    RSA_Data_Enc.TrialInfo{sub} = tmp_strct_enc.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;

    RSA_Data_Ret.Retrieval_Data{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.Retrieval_Data;
    RSA_Data_Ret.TrialInfo{sub} = tmp_strct_ret.(['RSA_Data_',Subj_names{sub}]).OCC.TrialInfo;

end
clear tmp_strct_enc tmp_strct_ret


% Plot Scattergrams of Real Data

ROI = {'OCC','TMP','FRT','CNT','PRT'};
r = 1;
N = 56;


figure('Pos', [326 71 1291 912])

for sbp = 1:9
    
    enc_samp = datasample(1:length(RSA_Data_Enc.TimeVec),1);
    ret_samp = datasample(1:length(RSA_Data_Ret.TimeVec),1);
    
    cur_data1 = nanmean(squeeze(RSA_Data_Enc.(ROI{r}).red16_Data(:,enc_samp,Perceptual_Mat_red16 > 0))',2);
    cur_data2 = nanmean(squeeze(RSA_Data_Ret.(ROI{r}).red16_Data(:,enc_samp,Perceptual_Mat_red16 > 0))',2);
    gr_WI = [reshape(cur_data1,[],1) reshape(cur_data2,[],1)];
    
    cur_data1 = nanmean(squeeze(RSA_Data_Enc.(ROI{r}).red16_Data(:,enc_samp,Perceptual_Mat_red16 < 0))',2);
    cur_data2 = nanmean(squeeze(RSA_Data_Ret.(ROI{r}).red16_Data(:,enc_samp,Perceptual_Mat_red16 < 0))',2);
    gr_BT = [reshape(cur_data1,[],1) reshape(cur_data2,[],1)];
    
    R_tr = [gr_WI; gr_BT];
    
    h = subplot(3,3,sbp);
    plot(gr_WI(:,1),gr_WI(:,2),'ro')
    hold on
    plot(gr_BT(:,1),gr_BT(:,2),'go')  
    hold on
    B = polyfit(R_tr(:,1), R_tr(:,2),1);
    ls1 = plot(linspace(min(R_tr(:,1)), max(R_tr(:,1)),N),B(2) + linspace(min(R_tr(:,1)), max(R_tr(:,1)),N)*B(1),'--m','linewidth',1.5);
    text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R_tr(:,1), R_tr(:,2))),'0.','.'),'Color','m')
    B = polyfit(R_tr(1:N,1), R_tr(1:N,2),1);
    ls2 = plot(linspace(min(R_tr(1:N,1)), max(R_tr(1:N,1)),N),B(2) + linspace(min(R_tr(1:N,1)), max(R_tr(1:N,1)),N)*B(1),'--r','linewidth',1.5);
    text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R_tr(1:N,1), R_tr(1:N,2))),'0.','.'),'Color','r')
    B = polyfit(R_tr(N+1:end,1), R_tr(N+1:end,2),1);
    ls3 = plot(linspace(min(R_tr(N+1:end,1)), max(R_tr(N+1:end,1)),N),B(2) + linspace(min(R_tr(N+1:end,1)), max(R_tr(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
    text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R_tr(N+1:end,1), R_tr(N+1:end,2))),'0.','.'),'Color','g')
    hold off
    [~,p_val] = ttest2(reshape(R_tr(1:N,:),N*2,1), reshape(R_tr(N+1:end,:),N*2,1));
    meas = (abs(norminv(p_val,0,1))/(atanh(corr(gr_WI(:,1), gr_WI(:,2))) - atanh(corr(gr_BT(:,1), gr_BT(:,2)))))*(atanh(corr(R_tr(:,1), R_tr(:,2))));
    title(sprintf('Enc: %1.3f - Ret: %1.3f - V_{corr} = %1.4f', RSA_Data_Enc.TimeVec(enc_samp), RSA_Data_Ret.TimeVec(ret_samp), meas),'fontsize',8)
    axis square
    set(h,'xlim',[0.4 0.65],'ylim',[0.4 0.65])
    axis square
end



N = 56;

Dur_Meth = zeros(2,1);

tic
X1 = kron(eye(23),[kron([1;0],ones(N,1)) kron([0;1],ones(N,1))]); % ones(length(y2),1)
for i = 1:10000
    
    dat1 = rand(112,23);
    dat2 = rand(112,23);
    y = dat1.*dat2;
    y2 = bsxfun(@times, kron(eye(23),ones(2*N,1)),y(:));
    SS_tot1 = sum(bsxfun(@minus, y, mean(y)).^2,1);%' * bsxfun(@minus, y, mean(y));
    %SS_tot1 = bsxfun(@minus, y, mean(y))' * bsxfun(@minus, y, mean(y));
    SS_res1 = sum((y2 - X1*((X1'*X1)\X1'*y2)).^2,1);
    %SS_res1 = (y2 - X1*((X1'*X1)\X1'*y2))' * (y2 - X1*((X1'*X1)\X1'*y2));
    meas1 = atanh(sqrt(1 - SS_res1./SS_tot1))';
    
end
Dur_Meth(1,1) = toc;

tic
X2 = [kron([1;0], ones(N,1)), kron([0;1], ones(N,1))];
for i = 1:10000
    
    dat1 = rand(112,23);
    dat2 = rand(112,23);
    meas2 = zeros(23,1);
    SS_tot2 = zeros(23,1);
    SS_res2 = zeros(23,1);
    for sub = 1:23
        y = tiedrank_(dat1(:,sub) .* dat2(:,sub),1);
        SS_tot2(sub) = sum((y - mean(y)).^2);
        SS_res2(sub) = sum((y - X2*((X2'*X2)\X2'*y)).^2);
        meas2(sub,1) = atanh(sqrt(1 - SS_res2(sub)/SS_tot2(sub)));
    end
    
end
Dur_Meth(2,1) = toc;





N = 56;

mu = [0 0; 3 3];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

figure('Pos', [326   500   619   483])
scatter(R(1:56,1),R(1:56,2),'ro','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
scatter(R(57:end,1),R(57:end,2),'go','MarkerFaceColor','g','MarkerEdgeColor','g')
set(gca,'xlim',[-3.5 6.5], 'ylim',[-3.5 6.5])
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--m','linewidth',1.5);
text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(R(:,1), R(:,2))),'0.','.'),'Color','m')
B = polyfit(R(1:N,1), R(1:N,2),1);
ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--r','linewidth',1.5);
text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(R(1:N,1), R(1:N,2))),'0.','.'),'Color','r')
B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(R(N+1:end,1), R(N+1:end,2))),'0.','.'),'Color','g')
hold off




