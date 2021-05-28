%% function [Cost] = ControlSimulator(U,P)
% Control simulator to evalute "true" cost of ido results. 
% IDO problem  "ControlSimulator.m" is solved with given optimal control 
% actions U, and optimal design and controller variables P.
% U(3,i)(energy from grid) that would incur negative SOC are forced to 
% purchase to maintain 0 SOC
%

function Cost1 = ControlSimulator(sol, P)
slv_time = tic;
X0 = [21;20;4;50]; %[inside temp, inside wall temp, wall core temp, SOC]

U0 = sol{3};
U = U0;

ah_kwh = 1000/12;

N = length(U);
nCluster = N/(24);

%% Periodicity
Q=24;
% Define time-varying constraints and disturbances
x_min = cell(Q,1);
x_max = cell(Q,1);
costs = cell(Q,1);

% grid pricing
peak_charge = .025;
offpeak = 0.01; %$/kWh
penalty = 1E-3;
Wxhc = ones(24,1)*0.1*penalty;

for k = 1:Q
    if k >= 8 && k <= 18
        x_min{k} = [21];
        x_max{k} = [26];   
        Wxhc(k) = [10E1]*penalty;
    else
        x_min{k} = [19];
        x_max{k} = [30]; 
    end

    %cost scheduale
    if k >8 && k<=20 % 8am:8pm = peak hrs
        costs{k} = peak_charge;   
    else
        costs{k} = offpeak;
    end 
end

%% simulate
b_size = P(2);
pv_size = P(1);

b_coef = 100/(b_size*ah_kwh);

%CAPEX params
pv_cost = 440; % $/m^2
b_cost = .136 *(12); %$/Ah

U(3,:) = sol{3};
X(1,:) = sol{2};

Cost = (pv_cost*pv_size + b_size*ah_kwh*b_cost)/(10); % CAPEX amortized monthly;
CAPEX = Cost;
if nCluster == 5
    for i = 1:nCluster
        k = (i-1)*24+1;
        j = i*24;
        Cost = Cost + U(3,k:j)*([costs{1:Q}])'*ClustProb(i); %energy purchased/sold to grid
        vhc = (max(X(1,k+1:j+1)-[x_max{:}],0))';
        vlc = (max([x_min{:}]-X(1,k+1:j+1),0))'; 
        Cost = Cost + [Wxhc(1:Q)]'*vhc*ClustProb(i); %high temp violations      
        Cost = Cost + [Wxhc(1:Q)]'*vlc*ClustProb(i); %low temp violations
%         if sum(vhc>5) >=1 %additional penalties for severe violations
% %             Cost = Cost+1E4*sum(vhc>5)*ClustProb(i) ;
%         elseif  sum(vlc>5) >=1
% %             Cost = Cost+1E4*sum(vlc>5)*ClustProb(i) ;
%         end

    end
    
elseif nCluster == 365
    Cost = Cost + U(3,:)*repmat([costs{1:Q}],1,nCluster)'; %energy purchased/sold to grid
    Ecost = Cost - CAPEX;
    vhc = (max(X(1,2:end)-repmat([x_max{mod(1:Q,Q)+1}],1,nCluster),0))'; %5 violation magnitude
    vlc = (max(repmat([x_min{mod(1:Q,Q)+1}],1,nCluster)-X(1,2:end),0))';
    Cost = Cost + repmat([Wxhc(1:Q)]',1,nCluster)*vhc; %high temp violations      
    Cost = Cost + repmat([Wxhc(1:Q)]',1,nCluster)*vlc; %low temp violations
%     if sum(vhc>5) >=1 %additional penalties for severe violations
% %         Cost = Cost+1E4*sum(vlc>5) ;
%     elseif  sum(vlc>5) >=1
% %         Cost = Cost+1E4*sum(vlc>5);
%     end
    Vcost = Cost - CAPEX - Ecost;
end
Cost1 = -Cost(end)/1E4; % scale objective

if 0 %save run data
    
    if 0% Case Study 1
        x1name = '../Output_Data/IDO_x1.mat';
        x2name = '../Output_Data/IDO_x2.mat';
        x3name = '../Output_Data/IDO_x3.mat';
        x4name = '../Output_Data/IDO_x4.mat';
        u1name = '../Output_Data/IDO_u1.mat';
        u2name = '../Output_Data/IDO_u2.mat';
        u3name = '../Output_Data/IDO_u3.mat';
        saveData(X(1,:), x1name, 'x1')
        saveData(X(2,:), x2name, 'x2')
        saveData(X(2,:), x3name, 'x3')
        saveData(X(4,:), x4name, 'x4')
        saveData(U(1,:), u1name, 'u1')
        saveData(U(2,:), u2name, 'u2')
        saveData(U(3,:), u3name, 'u3')     
        saveData(slv_time, 'IDO_solveTime', 'IDO_slv_times')
    elseif 1 %Case Study 2
        x1name = '../Output_Data/CLV_x1_1h_(1).mat';
        x2name = '../Output_Data/CLV_x2_1h_(1).mat';
        x3name = '../Output_Data/CLV_x3_1h_(1).mat';
        x4name = '../Output_Data/CLV_x4_1h_(1).mat';
        u1name = '../Output_Data/CLV_u1_1h_(1).mat';
        u2name = '../Output_Data/CLV_u2_1h_(1).mat';
        u3name = '../Output_Data/CLV_u3_1h_(1).mat';
        saveData(X(1,:), x1name, 'x1_1h_1')
        saveData(X(2,:), x2name, 'x2_1h_1')
        saveData(X(2,:), x3name, 'x3_1h_1')
        saveData(X(4,:), x4name, 'x4_1h_1')
        saveData(U(1,:), u1name, 'u1_1h_1')
        saveData(U(2,:), u2name, 'u2_1h_1')
        saveData(U(3,:), u3name, 'u3_1h_1')
        saveData(slv_time, 'CLV_solveTime', 'CLV_slv_times')        
    end
end

%% plot
if 0
    LB = [repmat([x_min{:}],1,nCluster);ones(1,Q*nCluster)*10];
    UB = [repmat([x_max{:}],1,nCluster);ones(1,Q*nCluster)*95];
    figure(); hold on;
    title(P)
    for i = 1:2
        k = [1,4];
        subplot(2,1,i); hold on;
        plot(0:(Q)*nCluster,X(k(i),:),'-or','linewidth',2)        
        plot(1:(Q)*nCluster,LB(i,:),'k--','linewidth',2)
        plot(1:(Q)*nCluster,UB(i,:),'k--','linewidth',2)
        set(gcf,'color','w');
        set(gca,'FontSize',16)
        ylabel("X_"+i)
        xlabel('time [hrs]')
    end

    yl = [0 1500;0 1500;-1500 1500];
    figure(); hold on;
    title(P)
    for i = 1:3
        subplot(3,1,i); hold on;
        plot(1:Q*nCluster,U(i,:),'-or','linewidth',2)
        set(gcf,'color','w');
        set(gca,'FontSize',16)
        ylim(yl(i,:))
        ylabel("u_"+i)
        xlabel('time [hrs]')
        box on
    end
end


end