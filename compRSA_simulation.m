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



mus1 = -2:1:2;
mus2 = -2:1:2;

mu = [0 0; 0 0];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
NMeas_Cell = cell(length(mus1),length(mus2));
for m1 = 1:length(mus1)
    for m2 = 1:length(mus2) 

        corrs1 = -0.9:0.1:0.9;
        corrs2 = -0.9:0.1:0.9;

        NMeas_Cell{m1,m2} = zeros(length(corrs1),length(corrs2));
        for c1 = 1:length(corrs1)
            for c2 = 1:length(corrs2)
                mu(2,:) = [mus1(m1) mus2(m2)];
                sigma(1,2,1) = corrs1(c1); sigma(2,1,1) = corrs1(c1);
                sigma(1,2,2) = corrs2(c2); sigma(2,1,2) = corrs2(c2);
                R = [];
                R = [R; mvnrnd(mu(1,:),sigma(:,:,1),0.5*N*(0.5*N -1))];
                R = [R; mvnrnd(mu(2,:),sigma(:,:,2),0.5*N*(0.5*N -1))];
                
                cur_data = tiedrank_(R,1);
                cor_all = fast_corr(cur_data(:,1), cur_data(:,2));
                cur_data = tiedrank_(R(1:56,:),1);
                cor_wi  = fast_corr(cur_data(:,1), cur_data(:,2));
                cur_data = tiedrank_(R(57:112,:),1);
                cor_bt  = fast_corr(cur_data(:,1), cur_data(:,2));
                NMeas_Cell{m1,m2}(c1,c2) = (atanh(cor_wi) - atanh(cor_bt))*atanh(cor_all);
            end
        end
        
    end
end


figure('Pos',[164 67 1504 926])
Data = NMeas_Cell';
tmp_max = cell2mat(cellfun(@(x) max(x(:)), Data,'UniformOutput',0));
tmp_min = cell2mat(cellfun(@(x) min(x(:)), Data,'UniformOutput',0));
c_limit = [-0.6 0.6];%[min(tmp_min(:))+0.9 max(tmp_max(:))-1];
phase_shift = [1.2 0.6 0 -0.6 -1.2];
for sbp = 1:25
    subplot(5,5,sbp)
    imagesc(Data{sbp})
    set(gca,'clim',c_limit, 'xtick',[], 'ytick',[])
    axis square
    axes('pos',[0.154 0.774 0.0772 0.025],'xlim',[-15 155],'ylim',[-15 15],'visible','off')
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
    axes('pos',[0.135 0.8 0.021 0.127],'xlim',[-15 15],'ylim',[-15 155],'visible','off')
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
end 


phase_shift = [1.2 0.6 0 -0.6 -1.2];
figure
hold on
for sbp = 1:5
    t = linspace(0,2*pi) ;
    a = 15 ; b = 15 ;
    x = a*cos(t+phase_shift(sbp)) ;
    y = b*sin(t) ;
    plot(x+35*(sbp-1),y,'k','linewidth',3)
end
axis equal
box off
set(gca,'xlim',[-15 155],'ylim',[-15 15],'visible','off')