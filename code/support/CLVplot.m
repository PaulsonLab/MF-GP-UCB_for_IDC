%% function [] = CLVplot(choice)
% plot trajectories for closed-loop verification
% Choice 0: get table values for Mean cost and CI of 4 cases
% Choice 1: plot 1 year long trajectories
% Choice 2: plot nd day trajectory from nd winter and summer days
% Choice 3: scatter nd day trajectory from nd winter and summer days 
% Choice 4: heatmap of trajectory (hour vs day)
% Choice 5: scatter of control trajectory
% Choice 6: heatmap of control trajectory
%

function [] = CLVplot(choice)
%% load states trajectories
addpath('../Output_Data')
%load X1(Temp) trajectory - best values (lowest cost)
load('CLV_x1_1h_(1).mat')       % low var solution on low var noise
load('CLV_x1_2h_(1).mat')       % low var solution on high var noise
load('CLV_x1_1h_(2).mat')       % high var solution on low var noise
load('CLV_x1_2h_(2).mat')       % high var solution on High var noise

%load X4(SOC) trajectory - best values (lowest cost)
load('CLV_x4_1h_(1).mat')       % low var solution on low var noise
load('CLV_x4_2h_(1).mat')       % low var solution on high var noise
load('CLV_x4_1h_(2).mat')       % high var solution on low var noise
load('CLV_x4_2h_(2).mat')       % high var solution on High var noise

%% load inputs trajectories
%load U1(heating) 
load('CLV_u1_1h_(1).mat')       % low var solution on low var noise
load('CLV_u1_2h_(1).mat')       % low var solution on high var noise
load('CLV_u1_1h_(2).mat')       % high var solution on low var noise
load('CLV_u1_2h_(2).mat')       % high var solution on High var noise

%load U2(cooling) 
load('CLV_u2_1h_(1).mat')       % low var solution on low var noise
load('CLV_u2_2h_(1).mat')       % low var solution on high var noise
load('CLV_u2_1h_(2).mat')   % high var solution on low var noise
load('CLV_u2_2h_(2).mat')   % high var solution on High var noise

%load U3(grid) 
load('CLV_u3_1h_(1).mat')       % low var solution on low var noise
load('CLV_u3_2h_(1).mat')       % low var solution on high var noise
load('CLV_u3_1h_(2).mat')   % high var solution on low var noise
load('CLV_u3_2h_(2).mat')   % high var solution on High var noise

%% compile solutions for easy plotting
% compile state trajectories
X1 = {x1_1h_1 x1_2h_1(511-270+1:end,:) x1_1h_2 x1_2h_2 };
X4 = {x4_1h_1 x4_2h_1(511-270+1:end,:) x4_1h_2 x4_2h_2 };
% compile control trajectories
U1 = {u1_1h_1 u1_2h_1(511-270+1:end,:) u1_1h_2 u1_2h_2 };
U2 = {u2_1h_1 u2_2h_1(511-270+1:end,:) u2_1h_2 u2_2h_2 };
U3 = {u3_1h_1 u3_2h_1(511-270+1:end,:) u3_1h_2 u3_2h_2 };

%% plotting set-up
Legend_list = {'f_1(w_1)' 'f_1(w_2)' 'f_2(w_1)' 'f_2(w_2)' };
t = linspace(0,365,8761);
nd = 3;

xhcLb = [19*ones(1,8),21*ones(1,10),19*ones(1,6)]; 
xhcUb = [30*ones(1,8),26*ones(1,10),30*ones(1,6)];

%% compute Costs
if choice ==0
    load('dh1.mat')
    load('dh2.mat')

    Cost = zeros(270,4);
    for i = 1:270
        for j = 1:4
            sol = {X1{j}(i,:),X4{j}(i,:),U3{j}(i,:)};
            if i<3
                d = dh1;
            else
                d=dh2;
            end
            Cost(i,j) = ControlSimulator2(sol,d);
        end
    end
    meanCost = -mean(Cost);
    worstCost = -min(Cost);
    ciCost = std(Cost)/sqrt(270);

    fprintf('Mean cost +/- CI, and worst cost for %2i Bayesian runs: \n',270)
    fprintf('f_1(w_1) %2f +/- %2f;   %2f\n',meanCost(1),ciCost(1),worstCost(1))
    fprintf('f_1(w_2) %2f +/- %2f;   %2f\n',meanCost(2),ciCost(2),worstCost(2))
    fprintf('f_2(w_1) %2f +/- %2f;   %2f\n',meanCost(3),ciCost(3),worstCost(3))
    fprintf('f_2(w_2) %2f +/- %2f;   %2f\n',meanCost(4),ciCost(4),worstCost(4))
    
    
    % compute violations
    pk_hrs = zeros(24,365);
    pk_hrs(8:18,:) = 1;
    pk_hrs = reshape(pk_hrs, 1, []);
    pk_idx =find(pk_hrs==1);

    x1 = X1{1}(:,pk_idx+1);
    x1_= mean(x1);
    x2 = X1{2}(:,pk_idx+1);
    x2_ = mean(x2);
    x3 = X1{3}(:,pk_idx+1);
    x3_ = mean(x3); 
    x4 = X1{4}(:,pk_idx+1);
    x4_ = mean(x4);
    
    
    Viol1 = (x1>26) + (x1<21);
    Viol2 = (x2>26) + (x2<21);
    Viol3 = (x3>26) + (x3<21);
    Viol4 = (x4>26) + (x4<21); 
    
    Viol1_total = sum(sum(Viol1));
    Viol2_total = sum(sum(Viol2));
    Viol3_total = sum(sum(Viol3));
    Viol4_total = sum(sum(Viol4)); 
    
    Viol1_magnitudes_avg = sum(sum(max(x1-26,0)+abs(min((x1-21),0))))/Viol1_total;
    Viol2_magnitudes_avg = sum(sum(max(x2-26,0)+abs(min((x2-21),0))))/Viol2_total;
    Viol3_magnitudes_avg = sum(sum(max(x3-26,0)+abs(min((x3-21),0))))/Viol3_total;
    Viol4_magnitudes_avg = sum(sum(max(x4-26,0)+abs(min((x4-21),0))))/Viol4_total;
    
    fprintf('Mean frequency of peak violations +/- CI, and max violation for %2i Bayesian runs: \n',270)
    fprintf('f_1(w_1) %2f +/- %2f;   %2f\n',...
        Viol1_total/(length(pk_idx)*270),std(std(Viol1))/(270^.5),Viol1_magnitudes_avg)
    fprintf('f_1(w_2) %2f +/- %2f;   %2f\n',...
        Viol2_total/(length(pk_idx)*270),std(std(Viol2))/(270^.5),Viol2_magnitudes_avg)
    fprintf('f_2(w_1) %2f +/- %2f;   %2f\n',...
        Viol3_total/(length(pk_idx)*270),std(std(Viol3))/(270^.5),Viol3_magnitudes_avg)
    fprintf('f_2(w_2) %2f +/- %2f;   %2f\n',...
        Viol4_total/(length(pk_idx)*270),std(std(Viol4))/(270^.5),Viol4_magnitudes_avg)
    
 

end
        
%% Choice 1
if choice == 1

    
    figure(1); 
    for i = 1:4
        subplot(4,1,i); hold on;
        plot(t, mean(X1{i}), 'b')
        plot(t, X1{i}, ':b')
        plot(t(2:end), repmat(xhcLb,1,365), ':k', 'Linewidth', 2)
        plot(t(2:end), repmat(xhcUb,1,365), ':k', 'Linewidth', 2)
        title(Legend_list{i})    
        xlabel('Day')
        ylabel('Temp ^o C')
    end
    sgtitle('X_1 Trajectory')

    figure(2);
    for i = 1:4
        subplot(4,1,i); hold on;
        plot(t, mean(X4{i}),'b')
        plot(t, (X4{i}),':b')
        plot(t, zeros(1,24*365+1), ':k', 'Linewidth', 2)
        plot(t, ones(1,24*365+1)*100, ':k', 'Linewidth', 2)
        title(Legend_list{i})
        xlabel('Day')
        ylabel('SOC')
        ylim([0 100])
    end
    sgtitle('X_4 Trajectory')
end

%% Choice 2
if choice == 2    
	for i = 1:4
        figure();  hold on;
        for j = 1:270
            p1 = plot(t(1:24*nd), X1{i}(j,5*24:24*(nd+5)-1),':b');
            p2 = plot(t(1:24*nd), X1{i}(j,180*24:24*(nd+180)-1),':r');
            p1.Color(4) = 0.2;
            p2.Color(4) = 0.2;
        end
        plot(t(1:24*nd), repmat(xhcLb,1,nd), ':k', 'Linewidth', 2)
        plot(t(1:24*nd), repmat(xhcUb,1,nd), ':k', 'Linewidth', 2)
        box on
%         xlabel('Day')
%         ylabel('Temp [^oC]')
        set(gcf,'color','w');
        set(gca,'FontSize',20)
        ylim([18 31])
        sname = "../../results3/CS2/X1"+Legend_list{i};
        saveas(gcf,sname,'pdf')
        saveas(gcf,sname,'fig')
    end

    for i = 1:4
        figure(); hold on;
        for j = 1:270
            p1 = plot(t(1:24*nd), X4{i}(j,5*24:24*(nd+5)-1),':b');
            p2 = plot(t(1:24*nd), X4{i}(j,180*24:24*(nd+180)-1),':r');        
            p1.Color(4) = 0.2;
            p2.Color(4) = 0.2;
        end
        plot(t(1:24*nd), zeros(1,24*nd), ':k', 'Linewidth', 2)
        plot(t(1:24*nd), ones(1,24*nd)*100, ':k', 'Linewidth', 2)
        box on
        set(gcf,'color','w');
        set(gca,'FontSize',20)
%         xlabel('Day')
%         ylabel('SOC [%]')
        ylim([0 100])
        sname = "../../results3/CS2/X4"+Legend_list{i};
        saveas(gcf,sname,'pdf')
        saveas(gcf,sname,'fig')
    end
    
%     figure(5); 
%     for i = 1:4
%         subplot(4,1,i); hold on;
%         p1 = plot(t(1:24*nd), X1{i}(:,362*24:24*(nd+362)-1),':b');
%         p2 = plot(t(1:24*nd), X1{i}(:,174*24:24*(nd+174)-1),':r');
%         plot(t(1:24*nd), repmat(xhcLb,1,nd), ':k', 'Linewidth', 2)
%         plot(t(1:24*nd), repmat(xhcUb,1,nd), ':k', 'Linewidth', 2)
%         title(Legend_list{i})    
%         xlabel('Day')
%         ylabel('Temp ^oC')
%     end
%     sgtitle('X_1 Trajectory(Bad days)')
% 
%     figure(6);
%     for i = 1:4
%         subplot(4,1,i); hold on;
%         p1 = plot(t(1:24*nd), X4{i}(:,362*24:24*(nd+362)-1),':b');
%         p2 = plot(t(1:24*nd), X4{i}(:,174*24:24*(nd+174)-1),':r');
%         plot(t(1:24*nd), zeros(1,24*nd), ':k', 'Linewidth', 2)
%         plot(t(1:24*nd), ones(1,24*nd)*100, ':k', 'Linewidth', 2)
%         title(Legend_list{i})
%         xlabel('Day')
%         ylabel('SOC')
%         ylim([0 100])
%     end
%     sgtitle('X_2 Trajectory(Bad days)')
    
end

%% Choice 3
if choice == 3
%     nd = 3;
%     figure(5); 
%     sgtitle('X_1 Trajectory')
%     for i = 1:8
%         subplot(4,2,i); hold on;
%         for j = 1:100
%             p1 = scatter(t(1:24*nd), X1{i}(j,1:24*nd),'MarkerFaceColor',...
%                 'b','MarkerEdgeColor','b');
%             p2 = scatter(t(1:24*nd), X1{i}(j,6*30*24+1:6*30*24+24*nd),...
%                 'MarkerFaceColor','r','MarkerEdgeColor','r');
%             p1.MarkerFaceAlpha = .1;	p1.MarkerEdgeAlpha = .1;
%             p2.MarkerFaceAlpha = .1;	p2.MarkerEdgeAlpha = .1;
%         end
%         plot(t(1:24*nd), repmat(xhcLb,1,nd), ':k', 'Linewidth', 2)
%         plot(t(1:24*nd), repmat(xhcUb,1,nd), ':k', 'Linewidth', 2)
%         title(Legend_list{i})  
%         ylim([10 40])
%         xlabel('Day')
%         ylabel('Temp ^o C')
%     end
% 
%     figure(6);
%     sgtitle('X_2 Trajectory')
%     for i = 1:8
%         subplot(4,2,i); hold on;
%         for j = 1:100
%             p1 = scatter(t(1:24*nd), X4{i}(j,1:24*nd),'MarkerFaceColor',...
%                 'b','MarkerEdgeColor','b');
%             p2 = scatter(t(1:24*nd), X4{i}(j,6*30*24+1:6*30*24+24*nd),...
%                 'MarkerFaceColor','r','MarkerEdgeColor','r');
%             p1.MarkerFaceAlpha = .1;	p1.MarkerEdgeAlpha = .1;
%             p2.MarkerFaceAlpha = .1;	p2.MarkerEdgeAlpha = .1;
%         end
%         plot(t(1:24*nd), zeros(1,24*nd), ':k', 'Linewidth', 2)
%         plot(t(1:24*nd), ones(1,24*nd)*100, ':k', 'Linewidth', 2)
%         title(Legend_list{i})
%         xlabel('Day')
%         ylabel('SOC')
%         ylim([0 100])
%     end
 end

%% Choice 4
if choice == 4
    figure(7)
    sgtitle('Mean Temperature Heat map')
    for i = 1:4
        subplot(4,1,i); 
        meanX = mean(X1{i});
        meanX(meanX(2:end)<repmat(xhcLb,1,365)) = nan;
        meanX(meanX(2:end)>repmat(xhcUb,1,365)) = nan;
        h = heatmap(reshape(meanX(2:end),24,365));    
        grid off
        colormap cool;   snapnow
        caxis([19, 30]);
        ylabel('Hour of day')
        xlabel('Day of year')
        title(Legend_list{i})
        %make special labels for heat map
        XLabels = 1:365;        YLabels = 1:24;
        % Convert each number in the array into a string
        CustomXLabels = string(XLabels);
        CustomYLabels = string(YLabels);
        % Replace all but the fifth elements by spaces
        CustomXLabels(mod(XLabels,30) ~= 0) = " ";
        CustomYLabels(mod(YLabels,4) ~= 0) = " ";
        % Set the 'XDisplayLabels' property of the heatmap 
        % object 'h' to the custom x-axis tick labels
        h.XDisplayLabels = CustomXLabels;
        h.YDisplayLabels = CustomYLabels;
    end

    figure(8)
    sgtitle('Mean SOC Heat map')
    for i = 1:4
        subplot(4,1,i); 
        meanX = mean(X4{i});
        meanX(meanX(2:end)<0) = nan;
        meanX(meanX(2:end)>100) = nan;
        h = heatmap(reshape(meanX(2:end),24,365));
        grid off
        colormap cool;  snapnow
        caxis([0, 100]);
        ylabel('Hour of day')
        xlabel('Day of year')
        title(Legend_list{i})
        %make special labels for heat map
        XLabels = 1:365;        YLabels = 1:24;
        % Convert each number in the array into a string
        CustomXLabels = string(XLabels);
        CustomYLabels = string(YLabels);
        % Replace all but the fifth elements by spaces
        CustomXLabels(mod(XLabels,30) ~= 0) = " ";
        CustomYLabels(mod(YLabels,4) ~= 0) = " ";
        % Set the 'XDisplayLabels' property of the heatmap 
        % object 'h' to the custom x-axis tick labels
        h.XDisplayLabels = CustomXLabels;
        h.YDisplayLabels = CustomYLabels;
    end
end

%% Choice 5
if choice == 5
    for i = 1:4
        figure();  hold on;
        for j = 1:270
            [a,b] = stairs(t(1:24*nd), U1{i}(j,5*24:24*(nd+5)-1)-U2{i}(j,5*24:24*(nd+5)-1),'b'); hold on;
            [c,d] = stairs(t(1:24*nd), U1{i}(j,180*24:24*(nd+180)-1)-U2{i}(j,180*24:24*(nd+180)-1),'r');hold on;
            p1 = plot(a,b,'b');
            p2 = plot(c,d,'r');
            p1.Color(4) = 0.2;
            p2.Color(4) = 0.2;
        end
        box on
        ylim([-1500 1500])
%         xlabel('Day')
%         ylabel('Net heating (U_1-U_2) [kWh]')
        set(gcf,'color','w');
        set(gca,'FontSize',20)
        sname = "../../results3/CS2/U1-U2_"+Legend_list{i};
        saveas(gcf,sname,'pdf')
        saveas(gcf,sname,'fig')
    end

    for i = 1:4
        figure(); hold on;
        for j = 1:270
            alpha_values = ones(24*nd,1)*.1;
            [a,b] = stairs(t(1:24*nd), U3{i}(j,5*24:24*(nd+5)-1),'b');
            [c,d] = stairs(t(1:24*nd), U3{i}(j,180*24:24*(nd+180)-1),'r');
            p1 = plot(a,b,'b');
            p2 = plot(c,d,'r');
            p1.Color(4) = 0.2;
            p2.Color(4) = 0.2;
        end
        box on
%         xlabel('Day')
%         ylabel('Energy from Grid (U_3) [kWh]')
        ylim([-1500 1500])
        set(gcf,'color','w');
        set(gca,'FontSize',20)
        sname = "../../results3/CS2/U3_"+Legend_list{i};
        saveas(gcf,sname,'pdf')
        saveas(gcf,sname,'fig')
    end
    
%     figure(11); 
%     for i = 1:4
%         subplot(4,1,i); hold on;
%         for j = 1:100
%             p1 = stairs(t(1:24*nd), U1{i}(j,362*24:24*(nd+362)-1)...
%                 -U2{i}(j,362*24:24*(nd+362)-1),'b');
%             p2 = stairs(t(1:24*nd), U1{i}(j,174*24:24*(nd+174)-1)...
%                 -U2{i}(j,177*24:24*(nd+177)-1),'r');
%         end
%         title(Legend_list{i})  
%         ylim([-1500 1500])
%         xlabel('Day')
%         ylabel('Net heating (U1-U2) ')
%     end
%     sgtitle('Net heating (U1-U2) (Bad days)')
% 
%     figure(12);
%     for i = 1:4
%         subplot(4,1,i); hold on;
%         for j = 1:100
%             p1 = stairs(t(1:24*nd), U3{i}(j,362*24:24*(nd+362)-1),'b');
%             p2 = stairs(t(1:24*nd), U3{i}(j,174*24:24*(nd+174)-1),'r');
%         end
%         title(Legend_list{i})
%         xlabel('Day')
%         ylabel('Energy from Grid')
%         ylim([-1500 1500])
%         sgtitle('Energy from Grid (U3) (Bad days)')
%     end
% end

%% Choice 6
if choice == 6
    figure(13)
    sgtitle('Mean Net Heating heatmap')
    for i = 1:4
        subplot(4,1,i); 
        meanX = mean(U1{i}-U2{i});
        h = heatmap(reshape(meanX,24,365));    
        grid off
        colormap cool;   snapnow
        caxis([-1500, 1500]);
        ylabel('Hour of day')
        xlabel('Day of year')
        title(Legend_list{i})
        %make special labels for heat map
        XLabels = 1:365;        YLabels = 1:24;
        % Convert each number in the array into a string
        CustomXLabels = string(XLabels);
        CustomYLabels = string(YLabels);
        % Replace all but the fifth elements by spaces
        CustomXLabels(mod(XLabels,30) ~= 0) = " ";
        CustomYLabels(mod(YLabels,4) ~= 0) = " ";
        % Set the 'XDisplayLabels' property of the heatmap 
        % object 'h' to the custom x-axis tick labels
        h.XDisplayLabels = CustomXLabels;
        h.YDisplayLabels = CustomYLabels;
    end

    figure(14)
    sgtitle('Mean energy from grid Heatmap')
    for i = 1:4
        subplot(4,1,i); 
        meanX = mean(U3{i});
        h = heatmap(reshape(meanX,24,365));
        grid off
        colormap cool;  snapnow
        caxis([-1500, 1500]);
        ylabel('Hour of day')
        xlabel('Day of year')
        title(Legend_list{i})
        %make special labels for heat map
        XLabels = 1:365;        YLabels = 1:24;
        % Convert each number in the array into a string
        CustomXLabels = string(XLabels);
        CustomYLabels = string(YLabels);
        % Replace all but the fifth elements by spaces
        CustomXLabels(mod(XLabels,30) ~= 0) = " ";
        CustomYLabels(mod(YLabels,4) ~= 0) = " ";
        % Set the 'XDisplayLabels' property of the heatmap 
        % object 'h' to the custom x-axis tick labels
        h.XDisplayLabels = CustomXLabels;
        h.YDisplayLabels = CustomYLabels;
    end
end



end