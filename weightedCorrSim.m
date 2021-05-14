%% Try out function callback

function weightedCorrSim(method)

mat_dim = 16;
n = 0.5*mat_dim*(0.5*mat_dim -1);
N = 2*n;

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
                Rnd_Data = [];
                Rnd_Data = [Rnd_Data; mvnrnd(mu(1,:),sigma(:,:,1),n)];
                Rnd_Data = [Rnd_Data; mvnrnd(mu(2,:),sigma(:,:,2),n)];
                
                Rnd_Data_tr = tiedrank_(Rnd_Data,1);
                cor_all = fast_corr(Rnd_Data(:,1), Rnd_Data(:,2));
                
                if(method == 1)
                    
                    %y = (Rnd_Data_tr(:,1) - mean(Rnd_Data_tr(:))) .* (Rnd_Data_tr(:,2) - mean(Rnd_Data_tr(:)));
                    y = tiedrank_(Rnd_Data(:,1) .* Rnd_Data(:,2),1);
                    X = [kron([1;0], ones(n,1)) kron([0;1], ones(n,1))];
                    SS_tot = (y - mean(y))' * (y - mean(y));
                    SS_res = (y - X*((X'*X)\X'*y))' * (y - X*((X'*X)\X'*y));
                    NMeas_Cell{m2,m1}(c2,c1) = atanh(sqrt(max(1 - SS_res/SS_tot,0))) * atanh(cor_all);
                
                else
                
                    X = [kron([1;0], ones(n,1)) kron([0;1], ones(n,1))];

                    group1 = Rnd_Data(logical(X(:,1)),:); group2 = Rnd_Data(logical(X(:,2)),:);
                    y1_m = mean(group1,1);
                    y2_m = mean(group2,1);
                    S1 = cov(group1);
                    S2 = cov(group2);
                    S = (((n-1)*S1 + (n-1)*S2)/(N-2))*((n*n)/(n + n));

                    T2 = (y1_m - y2_m)/inv(S)*(y1_m - y2_m)';
                    F = T2*((N - 3) / (2*(N - 2)));
                    %Rsq = sqrt(F/N)/(1+sqrt(F/N));
                    %R_OP = 1 - ((N - 3)/(N - 1 - 1))*(1 - Rsq)*(1 + (2*(1 - Rsq))/(N - 1 - 1));
                    %NMeas_Cell{m2,m1}(c2,c1) = atanh(sqrt(max(0,R_OP)))*atanh(cor_all);
                    NMeas_Cell{m2,m1}(c2,c1) = atanh(sqrt(sqrt(F/N)/(1+sqrt(F/N))))*atanh(cor_all);
                    %NMeas_Cell{m2,m1}(c2,c1) = f2z_bloc(F,1,N - 3)*atanh(cor_all);
                    
                end
                
            end
        end
        
    end
end

mus_long = [kron(mus2', ones(length(mus2),1)) kron(ones(length(mus1),1), mus1')];
cor_long = [kron(corrs1', ones(length(corrs2),1)) kron(ones(length(corrs1),1), corrs2')];


% Position in Pixels
figPos = [164 67 1504 926];
padd_bord_dat_x = 100;
padd_bord_y = 50;
min_input_area_size = 400;
data_area = [padd_bord_dat_x   padd_bord_y   figPos(3) - min_input_area_size - padd_bord_dat_x   figPos(4) - padd_bord_y];
padd_btw_x = 40;
padd_btw_y = 40;
data_plot_area = min(data_area([3 4])/5) - [padd_btw_x padd_btw_y];
data_plot_pos = zeros(25,2);
for i = 1:5
    for j = 1:5
        data_plot_pos((i-1)*5+j,1) = data_area(1) + (j-1)*(data_plot_area(1)+padd_btw_x);
        data_plot_pos((i-1)*5+j,2) = sum(data_area([2 4])) - (i)*(data_plot_area(2) + padd_btw_y);
    end
end
real_data_area = [padd_bord_dat_x   padd_bord_y   data_plot_pos(end,1)+(data_plot_area(1)+padd_btw_x)-padd_bord_dat_x  data_plot_pos(1,2)+(data_plot_area(2)+padd_btw_y)-padd_bord_y];
padd_bord_info_x = 120;
scatter_pos = [sum(real_data_area([1 3])) + padd_bord_info_x  0.5*(figPos(4) - (figPos(3) - sum(real_data_area([1 3])) - 2*padd_bord_info_x))  figPos(3) - sum(real_data_area([1 3])) - 2*padd_bord_info_x   figPos(3) - sum(real_data_area([1 3])) - 2*padd_bord_info_x];
hist_pos = [sum(real_data_area([1 3])) + padd_bord_info_x  real_data_area(2)+10  figPos(3) - sum(real_data_area([1 3])) - 2*padd_bord_info_x   0.6*(figPos(3) - sum(real_data_area([1 3])) - 2*padd_bord_info_x)];
colbar_size = [30 200];
colbar_pos = [sum(real_data_area([1 3])) + padd_bord_info_x + colbar_size(1)  figPos(4) - padd_bord_y - colbar_size(2)  colbar_size(1)  colbar_size(2)];
editmin_pos = [sum(real_data_area([1 3])) + padd_bord_info_x + colbar_size(1) + 100   colbar_pos(2) + colbar_pos(4)*0.75   60   20];
editmax_pos = [sum(real_data_area([1 3])) + padd_bord_info_x + colbar_size(1) + 180   colbar_pos(2) + colbar_pos(4)*0.75   60   20];
chngbtn_pos = [sum(real_data_area([1 3])) + padd_bord_info_x + colbar_size(1) + 110   editmax_pos(2) - 60   120   40];


f = figure('Visible','off','Position',figPos);
Data = NMeas_Cell';
c_limit = prctile(NMeas_Cell{3,3}(:),[0.1 99.9]); %[min(tmp_min(:))+0.9 max(tmp_max(:))-1];
phase_shift = [1.2 0.6 0 -0.6 -1.2];
for sbp = 1:25
    imgh(sbp) = axes('Units','pixels','pos',[data_plot_pos(sbp,:) data_plot_area],'visible','on','ZTickLabel',num2str(sbp));
    imagesc(Data{sbp});
    set(imgh(sbp),'clim',c_limit, 'xtick',[], 'ytick',[])
    axis square
    aspect = get(imgh(sbp),'PlotBoxAspectRatio');
    set(imgh(sbp),'Units','pixels');
    pos = get(imgh(sbp),'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(imgh(sbp),'Position',pos);
    
    cirxh(sbp) = axes('Units','pixels','pos',[pos(1) pos(2)-25 pos(3) 20],'xlim',[-15 155],'ylim',[-15 15],'visible','off');
    hold on
    for elp = 1:5
        t = linspace(0,2*pi) ;
        a = 15; b = 15;
        x = a*cos(t+phase_shift(elp)) ;
        y = b*sin(t) ;
        plot(x+40*(elp-1),y,'k','linewidth',2)
    end
    axis equal
    box off
    
    ciryh(sbp) = axes('Units','pixels','pos',[pos(1)-25 pos(2) 20 pos(4)],'xlim',[-15 15],'ylim',[-15 155],'visible','off');
    hold on
    for elp = 1:5
        t = linspace(0,2*pi) ;
        a = 15; b = 15;
        x = a*cos(t+phase_shift(elp)) ;
        y = b*sin(t) ;
        plot(x,y+40*(elp-1),'k','linewidth',2)
    end
    axis equal
    box off
    
    set(imgh(sbp),'Units','normalized')
    set(cirxh(sbp),'Units','normalized')
    set(ciryh(sbp),'Units','normalized')
    
end

%% Add Scatter Plot Axis
scatterh = axes('Units','pixels','pos',scatter_pos,'visible','on');
xlabel('X Dim'); ylabel('Y Dim');

%% Add Histogram Plot Axis
histh = axes('Units','pixels','pos',hist_pos,'visible','on');
ylabel('Count'); title(sprintf('M = %1.4f / SD = %1.4f',0,0))

%% Add Colorbar
colbarh = colorbar(imgh(13),'units','Pixels');
set(colbarh,'Position',colbar_pos)

%% Add Min and Max Clim Edit Box
uicontrol('Style','text','String','Min','Position',editmin_pos + [0 20 0 0]);
heditMin = uicontrol('Style','edit','String',c_limit(1),'Position',editmin_pos);
uicontrol('Style','text','String','Max','Position',editmax_pos + [0 20 0 0]);
heditMax = uicontrol('Style','edit','String',c_limit(2),'Position',editmax_pos);

%% Add Clim Change Button
hchngbtn = uicontrol('Style','pushbutton','String','Change c limits','Position',chngbtn_pos,'Callback',{@pushbutton_callback, heditMin, heditMax, imgh, colbarh});


%% Set Visibility of Function on
f.Visible = 'on';

%% Get an Input from Data

dcm_obj = datacursormode(f);
set(dcm_obj,'DisplayStyle','datatip',...
  'SnapToDataVertex','off','Enable','on','UpdateFcn',@(~,evt) chosenpoint_callback(evt))

%% Callback Functions
    function pushbutton_callback(~, ~, edit1, edit2, im, cb)

        min_clim = str2double(edit1.String);
        max_clim = str2double(edit2.String);
        for cc = 1:25
            set(im(cc),'clim',[min_clim max_clim])
        end
        set(cb,'ylim',[min_clim max_clim])

    end

    function [txt] = chosenpoint_callback(evt)
        
        cursorpos = get(evt, 'Position');
        for ax = 1:25
            if(evt.Target.Parent == imgh(ax))
                
                mu2 = [0 0; 0 0];
                sigma2 = cat(3,[1 0; 0 1], [1 0; 0 1]);
                
                mu2(2,:) = [mus_long(ax,2) mus_long(ax,1)];
                sigma2(1,2,1) = corrs1(cursorpos(1)); sigma2(2,1,1) = corrs1(cursorpos(1));
                sigma2(1,2,2) = corrs2(cursorpos(2)); sigma2(2,1,2) = corrs2(cursorpos(2));
                
                Rnd_Data2 = [];
                Rnd_Data2 = [Rnd_Data2; mvnrnd(mu2(1,:),sigma2(:,:,1),n)];
                Rnd_Data2 = [Rnd_Data2; mvnrnd(mu2(2,:),sigma2(:,:,2),n)];

                Rnd_Data2_tr = tiedrank_(Rnd_Data2,1);
                cor_all2 = fast_corr(Rnd_Data2(:,1), Rnd_Data2(:,2));
                
                axes(scatterh)
                s1 = plot(Rnd_Data2_tr(1:n,1),Rnd_Data2_tr(1:n,2),'o','Color',[1 0 0]);
                hold on
                s2 = plot(Rnd_Data2_tr(n+1:end,1),Rnd_Data2_tr(n+1:end,2),'o','Color',[0 1 0]);
                set(gca,'xlim',[min(Rnd_Data2_tr(:))-29 max(Rnd_Data2_tr(:))+22], 'ylim',[min(Rnd_Data2_tr(:))-29 max(Rnd_Data2_tr(:))+22]) %set(gca,'xlim',[-5 5], 'ylim',[-5 5], 'xtick', [],'ytick',[])
                B = polyfit(Rnd_Data2_tr(:,1), Rnd_Data2_tr(:,2),1);
                ls1 = plot(linspace(min(Rnd_Data2_tr(:,1)), max(Rnd_Data2_tr(:,1)),n),B(2) + linspace(min(Rnd_Data2_tr(:,1)), max(Rnd_Data2_tr(:,1)),n)*B(1),'--m','linewidth',1.5);
                t1 = text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(Rnd_Data2_tr(:,1), Rnd_Data2_tr(:,2))),'0.','.'),'Color','m');
                B = polyfit(Rnd_Data2_tr(1:n,1), Rnd_Data2_tr(1:n,2),1);
                ls2 = plot(linspace(min(Rnd_Data2_tr(1:n,1)), max(Rnd_Data2_tr(1:n,1)),n),B(2) + linspace(min(Rnd_Data2_tr(1:n,1)), max(Rnd_Data2_tr(1:n,1)),n)*B(1),'--r','linewidth',1.5);
                t2 = text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(Rnd_Data2_tr(1:n,1), Rnd_Data2_tr(1:n,2))),'0.','.'),'Color','r');
                if(((t2.Extent(1) >= t1.Extent(1) &  t2.Extent(1) <= sum(t1.Extent([1 3]))) | (sum(t2.Extent([1 3])) >= t1.Extent(1) &  sum(t2.Extent([1 3])) <= sum(t1.Extent([1 3])))) && ((t2.Extent(2) >= t1.Extent(2) &  t2.Extent(2) <= sum(t1.Extent([2 4])))))
                    t2.Position(2) = t1.Position(2) + t1.Extent(4);
                elseif(((t2.Extent(1) >= t1.Extent(1) &  t2.Extent(1) <= sum(t1.Extent([1 3]))) | (sum(t2.Extent([1 3])) >= t1.Extent(1) &  sum(t2.Extent([1 3])) <= sum(t1.Extent([1 3])))) && ((sum(t2.Extent([2 4])) >= t1.Extent(2) &  sum(t2.Extent([2 4])) <= sum(t1.Extent([2 4])))))
                    t2.Position(2) = t1.Position(2) - t2.Extent(4);
                end
                B = polyfit(Rnd_Data2_tr(n+1:end,1), Rnd_Data2_tr(n+1:end,2),1);
                ls3 = plot(linspace(min(Rnd_Data2_tr(n+1:end,1)), max(Rnd_Data2_tr(n+1:end,1)),n),B(2) + linspace(min(Rnd_Data2_tr(n+1:end,1)), max(Rnd_Data2_tr(n+1:end,1)),n)*B(1),'--g','linewidth',1.5);
                t3 = text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(Rnd_Data2_tr(n+1:end,1), Rnd_Data2_tr(n+1:end,2))),'0.','.'),'Color','g');
                if(((t3.Extent(1) >= t1.Extent(1) &  t3.Extent(1) <= sum(t1.Extent([1 3]))) | (sum(t3.Extent([1 3])) >= t1.Extent(1) &  sum(t3.Extent([1 3])) <= sum(t1.Extent([1 3])))) && ((t3.Extent(2) >= t1.Extent(2) &  t3.Extent(2) <= sum(t1.Extent([2 4])))))
                    t3.Position(2) = t1.Position(2) + t1.Extent(4);
                elseif(((t3.Extent(1) >= t1.Extent(1) &  t3.Extent(1) <= sum(t1.Extent([1 3]))) | (sum(t3.Extent([1 3])) >= t1.Extent(1) &  sum(t3.Extent([1 3])) <= sum(t1.Extent([1 3])))) && ((sum(t3.Extent([2 4])) >= t1.Extent(2) &  sum(t3.Extent([2 4])) <= sum(t1.Extent([2 4])))))
                    t3.Position(2) = t1.Position(2) - t3.Extent(4);
                end
                if(((t3.Extent(1) >= t2.Extent(1) &  t3.Extent(1) <= sum(t2.Extent([1 3]))) | (sum(t3.Extent([1 3])) >= t2.Extent(1) &  sum(t3.Extent([1 3])) <= sum(t2.Extent([1 3])))) && ((t3.Extent(2) >= t2.Extent(2) &  t3.Extent(2) <= sum(t2.Extent([2 4])))))
                    t3.Position(2) = t2.Position(2) + t2.Extent(4);
                elseif(((t3.Extent(1) >= t2.Extent(1) &  t3.Extent(1) <= sum(t2.Extent([1 3]))) | (sum(t3.Extent([1 3])) >= t2.Extent(1) &  sum(t3.Extent([1 3])) <= sum(t2.Extent([1 3])))) && ((sum(t3.Extent([2 4])) >= t2.Extent(2) &  sum(t3.Extent([2 4])) <= sum(t2.Extent([2 4])))))
                    t3.Position(2) = t2.Position(2) - t3.Extent(4);
                end
                hold off
                %pause(.6)
                set(s1.MarkerHandle,'FaceColorData',uint8(255*[1;0;0;0.2]),'EdgeColorData',uint8(255*[0;0;0;0]))
                set(s2.MarkerHandle,'FaceColorData',uint8(255*[0;1;0;0.2]),'EdgeColorData',uint8(255*[0;0;0;0]))
                
                
                hist_data = zeros(1000,1);
                for p = 1:1000
                    Rnd_Data2 = [];
                    Rnd_Data2 = [Rnd_Data2; mvnrnd(mu2(1,:),sigma2(:,:,1),n)];
                    Rnd_Data2 = [Rnd_Data2; mvnrnd(mu2(2,:),sigma2(:,:,2),n)];

                    Rnd_Data2_tr = tiedrank_(Rnd_Data2,1);
                    cor_all2 = fast_corr(Rnd_Data2(:,1), Rnd_Data2(:,2));

                    if(method == 1)

                        %y2 = (Rnd_Data2_tr(:,1) - mean(Rnd_Data2_tr(:))) .* (Rnd_Data2_tr(:,2) - mean(Rnd_Data2_tr(:)));
                        y2 = tiedrank_(Rnd_Data2(:,1) .* Rnd_Data2(:,2),1);
                        X2 = [kron([1;0], ones(n,1)) kron([0;1], ones(n,1))];
                        SS_tot2 = (y2 - mean(y2))' * (y2 - mean(y2));
                        SS_res2 = (y2 - X2*((X2'*X2)\X2'*y2))' * (y2 - X2*((X2'*X2)\X2'*y2));
                        hist_data(p) = atanh(sqrt(max(1 - SS_res2/SS_tot2,0))) * atanh(cor_all2);

                    else

                        X2 = [kron([1;0], ones(n,1)) kron([0;1], ones(n,1))];

                        group12 = Rnd_Data2(logical(X2(:,1)),:); group22 = Rnd_Data2(logical(X2(:,2)),:);
                        y1_m2 = mean(group12,1);
                        y2_m2 = mean(group22,1);
                        S12 = cov(group12);
                        S22 = cov(group22);
                        S_2 = (((n-1)*S12 + (n-1)*S22)/(N-2))*((n*n)/(n + n));

                        T22 = (y1_m2 - y2_m2)/inv(S_2)*(y1_m2 - y2_m2)';
                        F2 = T22*((N - 3) / (2*(N - 2)));
                        %Rsq2 = sqrt(F2/N)/(1+sqrt(F2/N));
                        %R_OP2 = 1 - ((N - 3)/(N - 1 - 1))*(1 - Rsq2)*(1 + (2*(1 - Rsq2))/(N - 1 - 1));
                        %hist_data(p) =  atanh(sqrt(max(0,R_OP2)))*atanh(cor_all2);
                        hist_data(p) =  atanh(sqrt(sqrt(F2/N)/(1+sqrt(F2/N))))*atanh(cor_all2);
                        %hist_data(p) =  f2z_bloc(F2,1,N - 3)*atanh(cor_all2);

                    end

                end
                
                axes(histh)
                hist(hist_data, 35)
                set(gca,'xlim',[-max(abs(hist_data)) max(abs(hist_data))].*1.5)
                title(sprintf('M = %1.4f / SD = %1.4f', mean(hist_data), std(hist_data)))
                
            end
        end
        txt = {num2str(cursorpos)};
        
    end

end




