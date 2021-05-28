function [J] = run_OC(decomp_by_cluster, d, Var,n_pred)
if n_pred == 24
    load('optimal_control_object.mat')
elseif n_pred == 12
	load('optimal_control_object_shortpredictions.mat')
else
    [Fctrl predictor] = build_OC(n_pred);
end
addpath(genpath('../../'))
load('humDist.mat')

X0 = [21;20;4;50]; %[inside temp, inside wall temp, wall core temp,SOC]
ah_kwh = 1000/12;

%integer variables
pv_size = d(1); % m^2; range = [0:540]
b_size = d(2) ;   %Ahr [0 15.5*12*1000*10]

beta_HC_u =  d(3); %[0 2*.4]
beta_HC_l =  d(4); %[0 2*.4]

if Var == 1
    var = [1,100,10];
elseif Var == 2
    var = [1,100,10].*2;

end

%% Draw Uncertainty
if decomp_by_cluster
    nCluster = 5;
    load("Clusterdata365_"+nCluster+".mat")
    Tod = (Centroids(:,1:24))';
    DHI = (Centroids(:,24+1:2*24))';
    DNI = (Centroids(:,2*24+1:end))';
    
    %collect cluster probabilities
    for i = 1:nCluster
       ClustProb(i) =   size(Cluster_Lump{i},2);
    end
    
    lag1 = 3; % days
    lag2 = lag1*24; %hours
else
    nCluster = 365;
    load('Uncertainty_Data.mat')
    end_day = 365*24*(4)-1 ; % 1 year sim ends on yr 4
    st_day = 365*24*(3);
    Tod = reshape(TempData(st_day:end_day,2),24,nCluster);
    DHI  = reshape(DHI_Data(st_day:end_day,2),24,nCluster);
    DNI  = reshape(DNI_Data(st_day:end_day,2),24,nCluster);
    
        lag1 = 3; % days
    lag2 = lag1*24; %hours
        %collect cluster probabilities
    for i = 1:nCluster
       ClustProb(i) =  1;% 30/nCluster;
    end
end

%% periodicity
Q=24;
% Define time-varying constraints and disturbances
x_min = cell(Q,1);
x_max = cell(Q,1);
costs = cell(Q,1);
OI = cell(Q,1); %Occucpancy Index

% grid pricing
peak_charge = .025;
offpeak = 0.01; %$/kWh
penalty = 1;
penalty2 = 1;
Wxhc = ones(24,1)*10^2*penalty;
Wxhc2 = ones(24,1)*10^2*penalty2;

for k = 1:Q
    if k >= 8 && k <= 18
        x_min{k} = [21]+ [beta_HC_l];
        x_max{k} = [26]- [beta_HC_u];   
        OI{k} = [1];
        Wxhc(k) = [10^3]*penalty;
        Wxhc2(k) = [10^3]*penalty2;

    else
        x_min{k} = [19]+ [beta_HC_l];
        x_max{k} = [30]- [beta_HC_u]; 
        OI{k}= [0];
    end

    %cost scheduale
    if k >8 && k<=20 % 8am:8pm = peak hrs
        costs{k} = peak_charge;   
    else
        costs{k} = offpeak;
    end 
end

b_coef = 100/(b_size*ah_kwh);

%% run simulation
Nsim = nCluster;
%n_pred = 24;
nu = 3;

X = zeros(4,Q,Nsim);
U = zeros(nu,Q,Nsim);

%CAPEX params
pv_cost = 440; % $/m^2
b_cost = .136 *(12); %$/Ah

J = (pv_cost*pv_size + b_size*b_cost*ah_kwh)/(10) ; % CAPEX amortized ;
capex = (pv_cost*pv_size + b_size*b_cost*ah_kwh)/(10);
profit =0; viol_cost = 0;

X(:,1,1) = X0;
hr=0;
tic
for n = 1:Nsim
    if n ~=1 % set Continuity between days
        X(:,1,n) = X(:,end,n-1);
    end

    Jd = 0; %daily cost

    for t = 1:Q
        
        % Update Measurment history & forcast  
        human_fcast = 30*([OI{mod(t:t+n_pred-1,Q)+1}]); % predict mean based on Occupancy Index
        Tod_fcast = [Tod(t:end,n);Tod(1:t-1,mod(n,nCluster)+1)] + randn(24,1)*var(1);
        DNI_fcast = max([DNI(t:end,n);DNI(1:t-1,mod(n,nCluster)+1)] + randn(24,1)*var(2),0);
        DHI_fcast = max([DHI(t:end,n);DHI(1:t-1,mod(n,nCluster)+1)] + randn(24,1)*var(3),0); 
        
        % evaluate controller
        [sol_curr, ef] = Fctrl{X(:,t,n), [x_min{mod(t:t+n_pred,Q)+1}],...
            [x_max{mod(t:t+n_pred,Q)+1}],[Wxhc(mod(t:t+n_pred-1,Q)+1)],[DNI_fcast(1:n_pred)]',[DHI_fcast(1:n_pred)]',...
            [human_fcast],[Tod_fcast(1:n_pred)]',[costs{mod(t:t+n_pred-1,Q)+1}],[d(1:2)]};
     
        U(:,t,n) = sol_curr(:,1);
        
        %actual disturbances
        Tod_k = Tod(t,n);
        dni_k = DNI(t,n);
        dhi_k = DHI(t,n);
        hum_k = humDist((n-1)*24 +t); % generate human disturbances

        % get next state 
        X(4,t+1,n) =X(4,t,n)+b_coef*(-(U(1,t,n)+ U(2,t,n)) + U(3,t,n) + dhi_k*pv_size*0.18)*ah_kwh;
        if X(4,t+1,n) <0 % set xb=0, solve for U(3,i); buys grid energy instead of negative charge
            U(3,t,n) =  ((1-X(4,t,n))/(b_coef*ah_kwh)+(U(1,t,n)+U(2,t,n))-pv_size*.18*dhi_k);
            X(4,t+1,n) = X(4,t,n) + b_coef*( -(U(1,t,n)+U(2,t,n))+U(3,t,n) +pv_size*.18*dhi_k)*ah_kwh;
        end
        if X(4,t+1,n)>100
           X(4,t+1,n)=100 ;
        end           
        X(1:3,t+1,n) = Building_dynamics(X(1:3,t,n),U(1:2,t,n), [Tod_k;dni_k;hum_k]);

        % calculate running cost
        Eps = zeros(2,1);
        Eps(1) = max(X([1],t+1,n)-[x_max{mod(t,Q)+1}]',0);
        Eps(2) = max([x_min{mod(t,Q)+1}]'-X([1],t+1,n),0);

        Jd = Jd + U(3,t,n)*costs{t}+ [Wxhc2(k),Wxhc2(k)]*Eps;
        profit = profit + U(3,t,n)*costs{t};
        viol_cost = viol_cost +[Wxhc2(k),Wxhc2(k)]*Eps;
        
        hr = hr+1;
    end
    J = J+Jd*(ClustProb(n));
end

%run controlSimulator for adjusted price
sol = {reshape(X(1:3,:,:),3,[]), reshape(X(4,:,:),1,[]), reshape(U,3,[])};
J = ControlSimulator(sol, d);

end