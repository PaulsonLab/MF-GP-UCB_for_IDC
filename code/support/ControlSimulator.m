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

%collect states and pop 25th hr of each day
xb1 = reshape(sol{2}',25, nCluster );
xb = [X0(4), reshape(xb1(2:end,:), 1,[])];
xt1 = sol{1};
xt2 = reshape(xt1(1,:),25, nCluster );
xt = [X0(1), reshape(xt2(2:end,:), 1,[])];

X_IDO = [xt;xb];

%% load disturbances
load('humDist.mat')
if nCluster == 5
    load("Clusterdata365_"+nCluster+".mat")
    Tod = reshape(Centroids(:,1:24)',1,[])';
    DHI = reshape(Centroids(:,24+1:2*24)',1,[])';
    DNI = reshape(Centroids(:,2*24+1:end)',1,[])';
    
    %collect cluster probabilities
    for i = 1:nCluster
       ClustProb(i) =   size(Cluster_Lump{i},2);
    end
    
    lag1 = 3; % days
    lag2 = lag1*24; %hours
elseif nCluster == 365
    load('Uncertainty_Data.mat')
    end_day = 365*24*(4)-1 ; % 1 year sim ends on yr 4
    st_day = 365*24*(3);
    Tod = (TempData(st_day:end_day,2));
    DHI  = (DHI_Data(st_day:end_day,2));
    DNI  = (DNI_Data(st_day:end_day,2));
    
    lag1 = 3; % days
    lag2 = lag1*24; %hours
        %collect cluster probabilities
    for i = 1:nCluster
       ClustProb(i) =  1;% 30/nCluster;
    end
end

W = [Tod,DNI,OI(1:24*nCluster)',DHI]';

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
X=zeros(4,N+1);
X(:,1) = X0;
for i = 1:N
    Xb = X(4,i) + b_coef*( -(U(1,i)+U(2,i))+U(3,i) +pv_size*.18*W(4,i))*ah_kwh;
    if Xb <0 % set xb=0, solve for U(3,i); buys grid energy instead of negative charge
        
        U(3,i) =  (1-X(4,i))/(b_coef*ah_kwh) + (U(1,i)+U(2,i))-pv_size*.18*W(4,i);
        Xb = X(4,i) + b_coef*( -(U(1,i)+U(2,i))+U(3,i) +pv_size*.18*W(4,i))*ah_kwh;
    end
    if Xb>100
       Xb=100 ;
    end
    
    Xs = Building_dynamics(X(1:3,i),U(1:2,i),W(1:3,i));
    
    X(:,i+1) = [Xs;Xb];

end

%CAPEX params
pv_cost = 440; % $/m^2
b_cost = .136 *(12); %$/Ah

% compute corrected cost
if 1 % compute cost with given values
    U = sol{3};
    X = X_IDO;
end

Cost = (pv_cost*pv_size + b_size*ah_kwh*b_cost)/(10); % CAPEX divided by expected years of opperation;
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

    end
    
elseif nCluster == 365
    Cost = Cost + U(3,:)*repmat([costs{1:Q}],1,nCluster)'; %energy purchased/sold to grid
    Ecost = Cost - CAPEX;
    vhc = (max(X(1,2:end)-repmat([x_max{mod(1:Q,Q)+1}],1,nCluster),0))'; %5 violation magnitude
    vlc = (max(repmat([x_min{mod(1:Q,Q)+1}],1,nCluster)-X(1,2:end),0))';
    Cost = Cost + repmat([Wxhc(1:Q)]',1,nCluster)*vhc; %high temp violations      
    Cost = Cost + repmat([Wxhc(1:Q)]',1,nCluster)*vlc; %low temp violations
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