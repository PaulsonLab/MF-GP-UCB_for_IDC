%% function [] = IDOplot(choice)
% plot trajectories for closed-loop verification 
% Choice 1: plot 1 year long trajectories
% Choice 2: plot nd day trajectory from nd winter and summer days
% Choice 3: scatter nd day trajectory from nd winter and summer days 
% Choice 4: heatmap of trajectory (hour vs day)
% Choice 5: scatter of control trajectory
% Choice 6: heatmap of control trajectory
%

function [] = IDOplot(choice)
%% load states trajectories
addpath('../Output_Data')
%load X1(Temp) trajectory - best values (lowest cost)

load('IDO_x1.mat')       

%load X4(SOC) trajectory - best values (lowest cost)
load('IDO_x4.mat')       

%% load inputs trajectories

%load U1(heating) 
load('IDO_u1.mat')      

%load U2(cooling) 
load('IDO_u2.mat')       

%load U3(grid) 
load('IDO_u3.mat')       
%% compile solutions for easy plotting
% compile state trajectories
X1 = {x1};
X4 = {x4};
% compile control trajectories
U1 = {u1};
U2 = {u2};
U3 = {u3};

%% plotting set-up
t = linspace(0,365,8761);

xhcLb = [19*ones(1,8),21*ones(1,10),19*ones(1,6)]; 
xhcUb = [30*ones(1,8),26*ones(1,10),30*ones(1,6)];

%% choice 0
if choice == 0 
    
    pk_hrs = zeros(24,365);
    pk_hrs(8:18,:) = 1;
    pk_hrs = reshape(pk_hrs, 1, []);
    pk_idx = find(pk_hrs==1);

    x1_pk = (x1(:,pk_idx+1));
    x1_mean_pk = mean(reshape(x1_pk, [],365));
    x1_mean_wk = mean(reshape(x1_mean_pk(2:end), [], 52)); 
        
    u4 = (u1-u2);
    u4_mean = mean(reshape(u4, [],365));
    u4_std = std(reshape(u4, [],365));
    u4_mean_wk = mean(reshape(u4_mean(:,2:end), [], 52)); 
    u4_std_wk = std(reshape(u4_mean(:,2:end), [], 52));
    u4_mean_mo = mean(reshape(u4_mean(:,6:end), [], 12));
    u4_std_mo = std(reshape(u4_mean(:,6:end), [], 12));
    
    u3_mean = mean(reshape(u3, [],365));
    u3_std = std(reshape(u3, [],365));
    u3_mean_wk = mean(reshape(u3_mean(:,2:end), [], 52));
    u3_std_wk = std(reshape(u3_mean(:,2:end), [], 52));    
    u3_mean_mo = mean(reshape(u3_mean(:,6:end), [], 12));
    u3_std_mo = std(reshape(u3_mean(:,6:end), [], 12));
    
%     figure();hold on;
%     stairs(u1-u2,'b', 'Linewidth', 2)
%     ylabel('net heating')
%     xlabel('hour')
%     set(gcf,'color','w');
%     set(gca,'FontSize',20)
%     box on
%    

    xpatch_s_day = [177, 180, 180, 177];
    xpatch_w_day = [5, 8, 8, 5];
    xpatch_s_wk = xpatch_s_day./7+1;
    xpatch_w_wk = xpatch_w_day./7+1;
    ypatch = [-800, -800, 200,200];

    
    figure(); hold on
    patch(xpatch_s_day,ypatch, 'r','EdgeColor','none','FaceAlpha',.3)
    patch(xpatch_w_day,ypatch, 'r','EdgeColor','none','FaceAlpha',.3)
    stairs(u4_mean,'b', 'Linewidth', 2)
    %errorbar(u4_mean, u4_std,'b', 'Linewidth', 1,'LineStyle','none')
    ylabel('average net heating')
    xlabel('day')  
    xlim([1 365])
    set(gcf,'color','w');
    set(gca,'FontSize',20)
    box on
    
    figure(); hold on
    patch(xpatch_s_wk,ypatch, 'r','EdgeColor','none','FaceAlpha',.3)
    patch(xpatch_w_wk,ypatch, 'r','EdgeColor','none','FaceAlpha',.3)    
    stairs(u4_mean_wk,'b', 'Linewidth', 2)
    %errorbar(u4_mean_wk, u4_std_wk,'b', 'Linewidth', 1,'LineStyle','none')
    ylabel('average net heating')
    xlabel('week')
    xlim([1 52])
    set(gcf,'color','w');
    set(gca,'FontSize',20)
    box on
    
%     figure(); hold on
%     stairs(u4_mean_mo,'b', 'Linewidth', 2)
%     errorbar(u4_mean_mo,u4_std_mo,'b', 'Linewidth', 1,'LineStyle','none')
%     ylabel('average net heating')
%     xlabel('month')
%     xlim([1 12])
%     set(gcf,'color','w');
%     set(gca,'FontSize',20)
%     box on
    
%     figure();hold on;
%     stairs(u3,'b', 'Linewidth', 2)
%     ylabel('Energy from grid')
%     xlabel('hour')
%     set(gcf,'color','w');
%     set(gca,'FontSize',20)
%     box on
%     
    figure(); hold on
    patch(xpatch_s_day,ypatch, 'r','EdgeColor','none','FaceAlpha',.3)
    patch(xpatch_w_day,ypatch, 'r','EdgeColor','none','FaceAlpha',.3)
    stairs(u3_mean,'b', 'Linewidth', 2)
    %errorbar(u3_mean,u3_std,'b', 'Linewidth', 1,'LineStyle','none' )
    ylabel('average energy from grid')
    xlabel('day') 
    xlim([1 365])
    set(gcf,'color','w');
    set(gca,'FontSize',20)
    box on
    
    figure(); hold on
    patch(xpatch_s_wk,ypatch, 'r','EdgeColor','none','FaceAlpha',.3)
    patch(xpatch_w_wk,ypatch, 'r','EdgeColor','none','FaceAlpha',.3)
    stairs(u3_mean_wk,'b', 'Linewidth', 2)
    %errorbar(u3_mean_wk, u3_std_wk,'b', 'Linewidth', 1,'LineStyle','none')
    ylabel('average energy from grid')
    xlabel('week')  
    xlim([1 52])
    set(gcf,'color','w');
    set(gca,'FontSize',20)
    box on
    
%     figure(); hold on
%     stairs(u3_mean_mo,'b', 'Linewidth', 2)
%     errorbar(u3_mean_mo,u3_std_mo,'b', 'Linewidth', 1,'LineStyle','none')
%     ylabel('average Energy from grid')
%     xlabel('month')  
%     xlim([1 12])
%     set(gcf,'color','w');
%     set(gca,'FontSize',20)
%     box on
    
end

%% Choice 1
if choice == 1
    figure(1); hold on;
    plot(t, X1{1}, 'b')
    plot(t(2:end), repmat(xhcLb,1,365), ':k', 'Linewidth', 2)
    plot(t(2:end), repmat(xhcUb,1,365), ':k', 'Linewidth', 2)
    xlabel('day')
    ylabel('temp ^o C')
    title('X_1 Trajectory')
    set(gcf,'color','w');
    set(gca,'FontSize',20)

    figure(2); hold on;
    plot(t, (X4{1}),'b')
    plot(t, zeros(1,24*365+1), ':k', 'Linewidth', 2)
    plot(t, ones(1,24*365+1)*100, ':k', 'Linewidth', 2)
    xlabel('day')
    ylabel('SOC')
    ylim([0 100])
    title('X_4 Trajectory')
    set(gcf,'color','w');
    set(gca,'FontSize',20)
end

%% Choice 2
if choice == 2    
    nd = 3;
    figure(3); hold on;
    p1 = plot(t(1:24*nd), X1{1}(:,24*177+2:24*180+1),'-r','Linewidth', 2);
    p2 = plot(t(1:24*nd), X1{1}(:,24*5+2:24*8+1),'-b','Linewidth', 2);
    plot(t(1:24*nd), repmat(xhcLb,1,nd), ':k', 'Linewidth', 2)
    plot(t(1:24*nd), repmat(xhcUb,1,nd), ':k', 'Linewidth', 2)
    xlabel('day')
    ylabel('temp ^oC')
    box on;
%     legend('summer','winter')
    set(gcf,'color','w');
    set(gca,'FontSize',20)

    figure(4); hold on;
    p1 = plot(t(1:24*nd), X4{1}(:,24*177+2:24*180+1),'-r','Linewidth', 2);
    p2 = plot(t(1:24*nd), X4{1}(:,24*5+2:24*8+1),'-b','Linewidth', 2);
    plot(t(1:24*nd), zeros(120,24*nd), ':k', 'Linewidth', 2)
    plot(t(1:24*nd), ones(120,24*nd)*100, ':k', 'Linewidth', 2)
    xlabel('day')
    ylabel('SOC')
    ylim([-5 105])
    box on;
%     legend('summer','winter')
    set(gcf,'color','w');
    set(gca,'FontSize',20)
end



%% Choice 4
if choice == 4
    
    figure(7);
    meanX = X1{1};
    meanX(meanX(2:end)<repmat(xhcLb,1,365)) = nan;
    meanX(meanX(2:end)>repmat(xhcUb,1,365)) = nan;
    h = heatmap(reshape(meanX(2:end),24,365));    
    grid off
    colormap cool;   snapnow
%     caxis([0, 50]);
    h.ColorLimits = [19, 30];
    ylabel('Hour of day')
    xlabel('Day of year')
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
    h.Title = 'Mean Temperature Heat map';
    set(gcf,'color','w');
    

    figure(8);
    meanX = X4{1};
    meanX(meanX(2:end)<0) = nan;
    meanX(meanX(2:end)>100) = nan;
    h = heatmap(reshape(meanX(2:end),24,365));
    grid off
    colormap cool;  snapnow
    caxis([0, 100]);
    ylabel('Hour of day')
    xlabel('Day of year')
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
    h.Title = 'Mean SOC Heat map';
    set(gcf,'color','w');
end

%% Choice 5
if choice == 5
    nd = 3;
    figure(9); hold on;
    p1 = stairs(t(1:24*nd), U1{1}(:,24*177+2:24*180+1)-U2{1}(:,24*177+2:24*180+1),'r','LineWidth',2);
    p2 = stairs(t(1:24*nd), U1{1}(:,24*5+2:24*8+1)...
        -U2{1}(:,24*5+2:24*8+1), 'b','LineWidth',2);
    ylim([-1550 1501])
    box on;
%     legend('summer','winter')
    xlabel('day')
    ylabel('net heating (U1-U2) [kWh]')
    set(gcf,'color','w');
    set(gca,'FontSize',20)

    figure(10); hold on;
    p1 = stairs(t(1:24*nd), U3{1}(:,24*177+2:24*180+1),'r','LineWidth',2);
    p2 = stairs(t(1:24*nd), U3{1}(:,24*5+2:24*8+1),'b','LineWidth',2);
%     legend('summer','winter')
    box on;
    xlabel('day')
    ylabel('energy from grid [kWh]')
    ylim([-1550 1501])
    set(gcf,'color','w');
    set(gca,'FontSize',20)
end

%% Choice 6
if choice == 6
    figure(11)
    meanX = (U1{1}-U2{1});
    h = heatmap(reshape(meanX,24,365));    
    grid off
    colormap cool;   snapnow
    caxis([-1500, 1500]);
    ylabel('Hour of day')
    xlabel('Day of year')
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
    h.Title = 'Net heating Heat map';
    set(gcf,'color','w');

    figure(12)
    meanX = (U3{1});
    h = heatmap(reshape(meanX,24,365));
    grid off
    colormap cool;  snapnow
    caxis([-1500, 1500]);
    ylabel('Hour of day')
    xlabel('Day of year')
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
    h.Title = 'Energy to Grid Heat map';
    set(gcf,'color','w');
end



end