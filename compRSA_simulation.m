%% Compare RSA Matrices Simulation

mat_dim = 16;
n = 0.5*mat_dim*(0.5*mat_dim -1);
N = 2*n;

method = 1;

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


weightedCorrSim(method)


%% Plot Correlations

mat_dim = 16;
n = 0.5*mat_dim*(0.5*mat_dim -1);
N = 2*n;

method = 1;

mus = [kron(ones(5,1),([-3 -1 0 1 3])') kron(([3 1 0 -1 -3])',ones(5,1))];

figure('Pos', [326 71 1291 912])

mu = [0 0; 0 0];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
for sbp = 1:25
    mu(2,:) = [mus(sbp,1) mus(sbp,2)];
    Rnd_Data = [];
    Rnd_Data = [Rnd_Data; mvnrnd(mu(1,:),sigma(:,:,1),n)];
    Rnd_Data = [Rnd_Data; mvnrnd(mu(2,:),sigma(:,:,2),n)];

    Rnd_Data_tr = tiedrank_(Rnd_Data,1);
    cor_all = fast_corr(Rnd_Data(:,1), Rnd_Data(:,2));
    
    h1 = subplot(5,5,sbp);
    s1 = plot(Rnd_Data_tr(1:n,1),Rnd_Data_tr(1:n,2),'o','Color',[1 0 0]);
    hold on
    s2 = plot(Rnd_Data_tr(n+1:end,1),Rnd_Data_tr(n+1:end,2),'o','Color',[0 1 0]);
    set(gca,'xlim',[min(Rnd_Data_tr(:))-29 max(Rnd_Data_tr(:))+22], 'ylim',[min(Rnd_Data_tr(:))-29 max(Rnd_Data_tr(:))+22], 'xtick', [],'ytick',[]) %set(gca,'xlim',[-5 5], 'ylim',[-5 5], 'xtick', [],'ytick',[])
    B = polyfit(Rnd_Data_tr(:,1), Rnd_Data_tr(:,2),1);
    ls1 = plot(linspace(min(Rnd_Data_tr(:,1)), max(Rnd_Data_tr(:,1)),n),B(2) + linspace(min(Rnd_Data_tr(:,1)), max(Rnd_Data_tr(:,1)),n)*B(1),'--m','linewidth',1.5);
    %annotation('textbox',[h1.Position(1)+h1.Position(3)*ls1.XData(end)/diff(get(gca,'xlim')) h1.Position(2)+h1.Position(4)*ls1.YData(end)/diff(get(gca,'ylim')) 0.06 0.03],'String',strrep(sprintf('r = %1.3f',corr(Rnd_Data_tr(:,1), Rnd_Data_tr(:,2))),'0.','.'),'Color','m','LineStyle','none')
    t1 = text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(Rnd_Data_tr(:,1), Rnd_Data_tr(:,2))),'0.','.'),'Color','m');
    B = polyfit(Rnd_Data_tr(1:n,1), Rnd_Data_tr(1:n,2),1);
    ls2 = plot(linspace(min(Rnd_Data_tr(1:n,1)), max(Rnd_Data_tr(1:n,1)),n),B(2) + linspace(min(Rnd_Data_tr(1:n,1)), max(Rnd_Data_tr(1:n,1)),n)*B(1),'--r','linewidth',1.5);
    %annotation('textbox',[h1.Position(1)+h1.Position(3)*ls2.XData(end)/diff(get(gca,'xlim')) h1.Position(2)+h1.Position(4)*ls2.YData(end)/diff(get(gca,'ylim')) 0.06 0.03],'String',strrep(sprintf('r = %1.3f',corr(Rnd_Data_tr(1:n,1), Rnd_Data_tr(1:n,2))),'0.','.'),'Color','r','LineStyle','none')
    t2 = text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(Rnd_Data_tr(1:n,1), Rnd_Data_tr(1:n,2))),'0.','.'),'Color','r');
    if(((t2.Extent(1) >= t1.Extent(1) &  t2.Extent(1) <= sum(t1.Extent([1 3]))) | (sum(t2.Extent([1 3])) >= t1.Extent(1) &  sum(t2.Extent([1 3])) <= sum(t1.Extent([1 3])))) && ((t2.Extent(2) >= t1.Extent(2) &  t2.Extent(2) <= sum(t1.Extent([2 4])))))
        t2.Position(2) = t1.Position(2) + t1.Extent(4);
    elseif(((t2.Extent(1) >= t1.Extent(1) &  t2.Extent(1) <= sum(t1.Extent([1 3]))) | (sum(t2.Extent([1 3])) >= t1.Extent(1) &  sum(t2.Extent([1 3])) <= sum(t1.Extent([1 3])))) && ((sum(t2.Extent([2 4])) >= t1.Extent(2) &  sum(t2.Extent([2 4])) <= sum(t1.Extent([2 4])))))
        t2.Position(2) = t1.Position(2) - t2.Extent(4);
    end
    B = polyfit(Rnd_Data_tr(n+1:end,1), Rnd_Data_tr(n+1:end,2),1);
    ls3 = plot(linspace(min(Rnd_Data_tr(n+1:end,1)), max(Rnd_Data_tr(n+1:end,1)),n),B(2) + linspace(min(Rnd_Data_tr(n+1:end,1)), max(Rnd_Data_tr(n+1:end,1)),n)*B(1),'--g','linewidth',1.5);
    %annotation('textbox',[h1.Position(1)+h1.Position(3)*ls3.XData(end)/diff(get(gca,'xlim')) h1.Position(2)+h1.Position(4)*ls3.YData(end)/diff(get(gca,'ylim')) 0.06 0.03],'String',strrep(sprintf('r = %1.3f',corr(Rnd_Data_tr(n+1:end,1), Rnd_Data_tr(n+1:end,2))),'0.','.'),'Color','g','LineStyle','none')
    t3 = text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(Rnd_Data_tr(n+1:end,1), Rnd_Data_tr(n+1:end,2))),'0.','.'),'Color','g');
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
    pause(.1)
    set(s1.MarkerHandle,'FaceColorData',uint8(255*[1;0;0;0.2]),'EdgeColorData',uint8(255*[0;0;0;0]))
    set(s2.MarkerHandle,'FaceColorData',uint8(255*[0;1;0;0.2]),'EdgeColorData',uint8(255*[0;0;0;0]))
    
    
    if(method == 1)

        %y = (Rnd_Data_tr(:,1) - mean(Rnd_Data_tr(:))) .* (Rnd_Data_tr(:,2) - mean(Rnd_Data_tr(:)));
        y = tiedrank_(Rnd_Data(:,1) .* Rnd_Data(:,2),1);
        X = [kron([1;0], ones(n,1)) kron([0;1], ones(n,1))];
        SS_tot = (y - mean(y))' * (y - mean(y));
        SS_res = (y - X*((X'*X)\X'*y))' * (y - X*((X'*X)\X'*y));
        meas = atanh(sqrt(max(1 - SS_res/SS_tot,0))) * atanh(cor_all);

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
        %meas =  atanh(sqrt(max(0,R_OP)))*atanh(cor_all);
        meas =  atanh(sqrt(sqrt(F/N)/(1+sqrt(F/N))))*atanh(cor_all);
        %meas =  f2z_bloc(F,1,N - 3)*atanh(cor_all);

    end
    title(sprintf('V_{corr} = %1.4f', meas),'fontsize',8)
    axis square
end




%% Plot Parameter Distributions

mat_dim = 16;
n = 0.5*mat_dim*(0.5*mat_dim -1);
N = 2*n;

method = 1;

n_perms = 1000;
mus1 = [-3 -1 0 1 3]; %-2:1:2;
mus2 = [3 1 0 -1 -3]; %2:-1:-2;

mu = [0 0; 0 0];
sigma = cat(3,[1 -.2; -.2 1], [1 -.4; -.4 1]);
NMeas_Cell = cell(length(mus2),length(mus1));
for m2 = 1:length(mus2)
    for m1 = 1:length(mus1) 
        
        mu(2,:) = [mus1(m1) mus2(m2)];
        NMeas_Cell{m2,m1} = zeros(n_perms,1);
        for p = 1:n_perms
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
                NMeas_Cell{m2,m1}(p) = atanh(sqrt(max(1 - SS_res/SS_tot,0))) * atanh(cor_all);

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
                %NMeas_Cell{m2,m1}(p) =  atanh(sqrt(max(0,R_OP)))*atanh(cor_all);
                NMeas_Cell{m2,m1}(p) =  atanh(sqrt(sqrt(F/N)/(1+sqrt(F/N))))*atanh(cor_all);
                %NMeas_Cell{m2,m1}(p) =  f2z_bloc(F,1,N - 3)*atanh(cor_all);

            end
            
        end
        
    end
end


figure('Pos', [326 71 1291 912])
Data = NMeas_Cell';
for sbp = 1:25
    
    h1 = subplot(5,5,sbp);
    hist(Data{sbp},35)
    set(gca,'xlim',[-2 2])
    title(sprintf('M = %1.4f / SD = %1.4f', mean(Data{sbp}), std(Data{sbp})),'fontsize',8)
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

TimeVec_Enc = RSA_Data_Enc.TimeVec;
TimeVec_Ret = RSA_Data_Ret.TimeVec;


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
    
    group1 = rand(112,23);
    group2 = rand(112,23);
    y = group1.*group2;
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
    
    group1 = rand(112,23);
    group2 = rand(112,23);
    meas2 = zeros(23,1);
    SS_tot2 = zeros(23,1);
    SS_res2 = zeros(23,1);
    for sub = 1:23
        y = tiedrank_(group1(:,sub) .* group2(:,sub),1);
        SS_tot2(sub) = sum((y - mean(y)).^2);
        SS_res2(sub) = sum((y - X2*((X2'*X2)\X2'*y)).^2);
        meas2(sub,1) = atanh(sqrt(1 - SS_res2(sub)/SS_tot2(sub)));
    end
    
end
Dur_Meth(2,1) = toc;





N = 56;

mu = [3 3; 0 0];
sigma = cat(3,[1 0; 0 1], [1 0; 0 1]);
%rng('default')  % For reproducibility
R = [];
R = [R; mvnrnd(mu(1,:),sigma(:,:,1),N)];
R = [R; mvnrnd(mu(2,:),sigma(:,:,2),N)];

figure('Pos', [326   200   619   483])
scatter(R(1:56,1),R(1:56,2),'ro','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
scatter(R(57:end,1),R(57:end,2),'go','MarkerFaceColor','b','MarkerEdgeColor','b')
set(gca,'xlim',[-3.5 6.5], 'ylim',[-3.5 6.5],'xtick','','ytick','')
B = polyfit(R(:,1), R(:,2),1);
ls1 = plot(linspace(min(R(:,1)), max(R(:,1)),N),B(2) + linspace(min(R(:,1)), max(R(:,1)),N)*B(1),'--k','linewidth',3.5);
%B = polyfit(R(1:N,1), R(1:N,2),1);
%ls2 = plot(linspace(min(R(1:N,1)), max(R(1:N,1)),N),B(2) + linspace(min(R(1:N,1)), max(R(1:N,1)),N)*B(1),'--b','linewidth',1.5);
%B = polyfit(R(N+1:end,1), R(N+1:end,2),1);
%ls3 = plot(linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N),B(2) + linspace(min(R(N+1:end,1)), max(R(N+1:end,1)),N)*B(1),'--r','linewidth',1.5);
hold off





%% Correlations of Real Data

measures = {'LDA','SVM','euclidian','euclidian_wcc'};
msr = 1;
r = 1;

tmp_struct = load(['Crosscompare_RSA_Matfiles/CrossComp_RSA_Enc_Ret_wiROI_noperm_msr1'],['CrossComp_RSA_',measures{msr}]);
Data = tmp_struct.(['CrossComp_RSA_',measures{msr}]).RSA_red16;
TimeX = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec1;
TimeY = tmp_struct.(['CrossComp_RSA_',measures{msr}]).TimeVec2;
ROI = tmp_struct.(['CrossComp_RSA_',measures{msr}]).ROI;

Corr_Data = squeeze(nanmean(Data.(ROI{r}).Corr,1));
Meth1_per_Data = squeeze(nanmean(Data.(ROI{r}).Meth1_per,1));
Meth1_sem_Data = squeeze(nanmean(Data.(ROI{r}).Meth1_sem,1));

interest_vals_corr = prctile(Corr_Data(:),[0.5 20 50 80 99.5]);
interest_vals_m1p = prctile(Meth1_per_Data(:),[0.5 20 50 80 99.5]);
interest_vals_m1s = prctile(Meth1_sem_Data(:),[0.5 20 50 80 99.5]);

min_vals_c = dsearchn(Corr_Data(:),interest_vals_corr');
min_vals_mp = dsearchn(Meth1_per_Data(:),interest_vals_m1p');
min_vals_ms = dsearchn(Meth1_sem_Data(:),interest_vals_m1s');
ScatterData = cell(5,3);
for i = 1:5
    [ind2, ind1] = find(Corr_Data == Corr_Data(min_vals_c(i)));
    ScatterData{i,1} = [squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Enc', TimeX(ind1)),Perceptual_Mat_red16(:) > 0),1))...
                        squeeze(nanmean(RSA_Data_Ret.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Ret', TimeY(ind2)),Perceptual_Mat_red16(:) > 0),1));...
                        squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Enc', TimeX(ind1)),Perceptual_Mat_red16(:) < 0),1))...
                        squeeze(nanmean(RSA_Data_Ret.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Ret', TimeY(ind2)),Perceptual_Mat_red16(:) < 0),1))];
    [ind2,ind1] = find(Meth1_per_Data == Meth1_per_Data(min_vals_mp(i)));
    ScatterData{i,2} = [squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Enc', TimeX(ind1)),Perceptual_Mat_red16(:) > 0),1))...
                        squeeze(nanmean(RSA_Data_Ret.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Ret', TimeY(ind2)),Perceptual_Mat_red16(:) > 0),1));...
                        squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Enc', TimeX(ind1)),Perceptual_Mat_red16(:) < 0),1))...
                        squeeze(nanmean(RSA_Data_Ret.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Ret', TimeY(ind2)),Perceptual_Mat_red16(:) < 0),1))];
    [ind2,ind1] = find(Meth1_sem_Data == Meth1_sem_Data(min_vals_ms(i)));
    ScatterData{i,3} = [squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Enc', TimeX(ind1)),Perceptual_Mat_red16(:) > 0),1))...
                        squeeze(nanmean(RSA_Data_Ret.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Ret', TimeY(ind2)),Perceptual_Mat_red16(:) > 0),1));...
                        squeeze(nanmean(RSA_Data_Enc.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Enc', TimeX(ind1)),Perceptual_Mat_red16(:) < 0),1))...
                        squeeze(nanmean(RSA_Data_Ret.(ROI{r}).red16_Data(:,dsearchn(TimeVec_Ret', TimeY(ind2)),Perceptual_Mat_red16(:) < 0),1))];
end

N = 56;

figure
ct = 1;
for sbp = 1:15
    curdat = ScatterData{ct};
    curdat_tr = tiedrank_(curdat,1);
    subplot(4,5,ct)
    s1 = plot(curdat_tr(1:N,1),curdat_tr(1:N,2),'o','Color',[1 0 0]);
    hold on
    s2 = plot(curdat_tr(N+1:end,1),curdat_tr(N+1:end,2),'o','Color',[0 1 0]);
    B = polyfit(curdat_tr(:,1), curdat_tr(:,2),1);
    ls1 = plot(linspace(min(curdat_tr(:,1)), max(curdat_tr(:,1)),N),B(2) + linspace(min(curdat_tr(:,1)), max(curdat_tr(:,1)),N)*B(1),'--m','linewidth',1.5);
    t1 = text(ls1.XData(end), ls1.YData(end),strrep(sprintf('r = %1.3f',corr(curdat_tr(:,1), curdat_tr(:,2))),'0.','.'),'Color','m');
    B = polyfit(curdat_tr(1:N,1), curdat_tr(1:N,2),1);
    ls2 = plot(linspace(min(curdat_tr(1:N,1)), max(curdat_tr(1:N,1)),N),B(2) + linspace(min(curdat_tr(1:N,1)), max(curdat_tr(1:N,1)),N)*B(1),'--r','linewidth',1.5);
    t2 = text(ls2.XData(end), ls2.YData(end),strrep(sprintf('r = %1.3f',corr(curdat_tr(1:N,1), curdat_tr(1:N,2))),'0.','.'),'Color','r');
    if(((t2.Extent(1) >= t1.Extent(1) &  t2.Extent(1) <= sum(t1.Extent([1 3]))) | (sum(t2.Extent([1 3])) >= t1.Extent(1) &  sum(t2.Extent([1 3])) <= sum(t1.Extent([1 3])))) && ((t2.Extent(2) >= t1.Extent(2) &  t2.Extent(2) <= sum(t1.Extent([2 4])))))
        t2.Position(2) = t1.Position(2) + t1.Extent(4);
    elseif(((t2.Extent(1) >= t1.Extent(1) &  t2.Extent(1) <= sum(t1.Extent([1 3]))) | (sum(t2.Extent([1 3])) >= t1.Extent(1) &  sum(t2.Extent([1 3])) <= sum(t1.Extent([1 3])))) && ((sum(t2.Extent([2 4])) >= t1.Extent(2) &  sum(t2.Extent([2 4])) <= sum(t1.Extent([2 4])))))
        t2.Position(2) = t1.Position(2) - t2.Extent(4);
    end
    B = polyfit(curdat_tr(N+1:end,1), curdat_tr(N+1:end,2),1);
    ls3 = plot(linspace(min(curdat_tr(N+1:end,1)), max(curdat_tr(N+1:end,1)),N),B(2) + linspace(min(curdat_tr(N+1:end,1)), max(curdat_tr(N+1:end,1)),N)*B(1),'--g','linewidth',1.5);
    t3 = text(ls3.XData(end), ls3.YData(end),strrep(sprintf('r = %1.3f',corr(curdat_tr(N+1:end,1), curdat_tr(N+1:end,2))),'0.','.'),'Color','g');
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
    axis equal
    pause(.1)
    set(s1.MarkerHandle,'FaceColorData',uint8(255*[1;0;0;0.5]),'EdgeColorData',uint8(255*[0;0;0;0]))
    set(s2.MarkerHandle,'FaceColorData',uint8(255*[0;1;0;0.5]),'EdgeColorData',uint8(255*[0;0;0;0]))
    set(gca,'xlim',[min(curdat_tr(:,1))-30 max(curdat_tr(:,1))+30],'ylim',[min(curdat_tr(:,2))-30 max(curdat_tr(:,2))+30])
    ct = ct + 1;
end
subplot(4,5,18)
hist(Corr_Data(:),30)





