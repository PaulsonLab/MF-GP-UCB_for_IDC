%% Integrated design and opperation for OCAM case study (no clusters)
% Build LP part of problem for CBO over NL parameters

function [controller_obj,controller_states] = build_IDO_yalmip(decomp_by_cluster, df_name, nCluster,sim_horz)
%% Set up values that dont change between days

addpath('../data')
nyears = 1; % nyears can be max 3 due to limits on current data set
st_day  = 365*24*3 ; %Start from year 3

if decomp_by_cluster
    load(df_name)
    nCluster = size(Centroids,1);
end

nu = 3;
ny = 1;

% State constraints 21-26C durring working hrs, 19-30 off hrs
xhcLb = [19*ones(1,9),21*ones(1,10),19*ones(1,6)]; %
xhcUb = [30*ones(1,9),26*ones(1,10),30*ones(1,6)]; %
xbLb = 10;
xbUb = 95;

% Input constraints
uLb = [0;0;-1500];
uUb =  1500*ones(1,nu)';

% Unit conversion
ah_kWh = 1000/12;

% grid pricing
peak_charge = .025;
offpeak = 0.01; %$/kWh

% Define time-varying constraints and disturbances
costs = zeros(24,1);
OI = zeros(24,1); 

penalty_test = 1;
Wxhc = ones(24,1)*10^2*penalty_test;

for k = 1:25
    %Occucpancy Index
    if k >= 8 && k <= 18 % 8am:6pm = occupancy effects 
        OI(k) = [30];
        Wxhc(k) = [10E3]*penalty_test;
    end
    %cost scheduale
    if k >8 && k<=20 % 8am:8pm = peak hrs
        costs(k) = peak_charge;   
    else
        costs(k) = offpeak;
    end
end

% initialze script / load data
if decomp_by_cluster
    %load cluster data
    end_day = st_day+nCluster*24;
    Tod = (Centroids(:,1:24))';
    DHI = (Centroids(:,24+1:2*24))';
    DNI = (Centroids(:,2*24+1:end))';
    
    %collect cluster probabilities
    for i = 1:nCluster
       ClustProb(i) =   size(Cluster_Lump{i},2);
    end
else
    load('Uncertainty_Data.mat')
    %end_day = 365*24*(3+nyears)-1 ; % 1 year sim ends on yr 4
    end_day = st_day + nCluster*24 -1;
    Tod = reshape(TempData(st_day:end_day,2),24,nCluster);
    DHI  = reshape(DHI_Data(st_day:end_day,2),24,nCluster);
    DNI  = reshape(DNI_Data(st_day:end_day,2),24,nCluster);
    
        %collect cluster probabilities
    for i = 1:nCluster
       ClustProb(i) =  1; % 30/nCluster;
    end
end

%% Build the IDO problem using YALMIP
tic
fprintf('Building IDO problem... ')
% define input and constraint vars
Umpc = sdpvar(repmat(nu,1,nCluster),repmat(24,1,nCluster),'full');
YHCmpc = sdpvar(repmat(3,1,nCluster),repmat(25,1,nCluster),'full');
YBmpc = sdpvar(repmat(1,1,nCluster),repmat(25,1,nCluster),'full');
obj = sdpvar(1,1);

%design vars
pv_size = sdpvar(1);
b_size  = sdpvar(1);

%Battery Model coefficient
b_coef = 100/(b_size*ah_kWh); 

%Building dynamics parameters
C1 = 9.356*10^5; C2 = 2.970*10^6; C3 = 6.695*10^5;
K1 = 16.48; K2 = 108.5; K3 = 5; K4 = 30.5; K5 = 23.04;
eta_c = 2;
eta_h = 4;

%CAPEX params
pv_cost = 440; % $/m^2
b_cost = .136 *(12); %$/Ah

%ammortization period\
terms = 10; %ammortize over ten years

objective = [];
objective = (pv_cost*pv_size + b_size*(ah_kWh)*b_cost)/(terms) ; % CAPEX
constraints = [];

%slack constraints
Eps = sdpvar(repmat(4,1,nCluster),repmat(24,1,nCluster),'full');
constraints = [constraints, 0<=[Eps{:}]];
%Eps = zeros(4,24,nCluster);

%build control part of problem
for j = 1:nCluster 
    if j ~=1  
        constraints = [constraints,YHCmpc{j}(:,1) == YHCmpc{j-1}(:,end)];
        constraints = [constraints,YBmpc{j}(1)  ==  YBmpc{j-1}(end)];
    end
    
    Jd = 0; %daily cost   
    for k = 1:24
        %input constraints
        constraints = [constraints, uLb <= Umpc{j}(:,k) <=uUb];

        %simplify variable names for readability
         wk = [Tod(k,j);DNI(k,j);OI(k)];
         
        %Temperature continuity constraints
        constraints = [constraints,YHCmpc{j}(1,k+1) == ...    
                3600*(1/C1)*((K1+K2)*(YHCmpc{j}(2,k)-YHCmpc{j}(1,k)) + K3*(wk(1)-YHCmpc{j}(1,k))+...
                K5*(YHCmpc{j}(3,k)-YHCmpc{j}(1,k))+wk(2)+Umpc{j}(1,k)*eta_h-Umpc{j}(2,k)*eta_c + wk(3))...
                + YHCmpc{j}(1,k)];
        constraints = [constraints,YHCmpc{j}(2,k+1) == ...
                3600*(1/C2)*((K1+K2)*(YHCmpc{j}(1,k)-YHCmpc{j}(2,k)) + wk(2)) + YHCmpc{j}(2,k)]; 
        constraints = [constraints,YHCmpc{j}(3,k+1) == ...
                3600*(1/C3)*(K5*(YHCmpc{j}(1,k)-YHCmpc{j}(3,k))+K4*(wk(1)-YHCmpc{j}(3,k))) + YHCmpc{j}(3,k)];     
           
        %SOC continuity constraints
        constraints = [constraints,YBmpc{j}(k+1) == ... %continuity
               YBmpc{j}(k)+b_coef*(-(Umpc{j}(1,k)+Umpc{j}(2,k)) +...
               Umpc{j}(3,k) + DHI(k,j)*pv_size*0.18)*ah_kWh];
           
        %SOC stat constraints   
        constraints = [constraints, (xbLb-Eps{j}(3,k) <= YBmpc{j}(k+1)...
                <= xbUb+Eps{j}(4,k)  )];
        %Temp state constraints  
        constraints = [constraints, (xhcLb(k+1)-Eps{j}(1,k) <= ...
               YHCmpc{j}(1,k+1) <= xhcUb(k+1)+Eps{j}(2,k) )];

        % Soft constraints
%         Eps(1,k,j) = max(YHCmpc{j}(1,k+1)-xhcUb(k+1),0);
%         Eps(2,k,j) = max(xhcLb(k+1)-YHCmpc{j}(1,k+1),0); 
%         
%         Eps(3,k,j) = max(YBmpc{j}(k+1)-xbUb,0);
%         Eps(4,k,j) = max(xbLb-YBmpc{j}(k+1),0); 
%         
%         constraints = [constraints, (xbLb-Eps(3,k,j) <= YBmpc{j}(k+1)...
%                 <= xbUb+Eps(4,k,j)  )];
%         constraints = [constraints, (xhcLb(k+1)-Eps(1,k,j) <= ...
%                YHCmpc{j}(1,k+1) <= xhcUb(k+1)+Eps(2,k,j) )];

        %Stage Cost
        %Jd = Jd + Umpc{j}(3,k)*costs(k)+ [Wxhc(k),Wxhc(k)]*Eps(1:2,k,j)+ Eps(3:4,k,j)'*diag([2E6,2E6])*Eps(3:4,k,j);
        Jd = Jd + Umpc{j}(3,k)*costs(k)+ [Wxhc(k),Wxhc(k)]*Eps{j}(1:2,k)+ Eps{j}(3:4,k)'*diag([2E6,2E6])*Eps{j}(3:4,k);
    end
    
  objective = objective+Jd*(ClustProb(j)); 
  %disp(j)

end
% Scale objective
objective = objective/1E4;

constraints = [constraints, obj == -objective];

%options = sdpsettings('solver','cplex','verbose',0,'debug',1);
options = sdpsettings('solver','gurobi','verbose',0,'debug',1) ; 
options.gurobi.OptimalityTol=1E-2; 
options.gurobi.FeasibilityTol = 1E-2;
options.gurobi.IntFeasTol=1E-2;

%The problem below evaluates IDO365 given I.Cs, disturbances, abd the
%battery size. Use blackbox optimization to solve the optimal design
controller_obj = optimizer(constraints, objective, options, ...
    {YHCmpc{1}(:,1),YBmpc{1}(1),[pv_size,b_size]},{obj});
controller_states = optimizer(constraints, objective, options, ...
    {YHCmpc{1}(:,1),YBmpc{1}(1),[pv_size,b_size]},{[YHCmpc{:}],[YBmpc{:}],[Umpc{:}]});
controller_eps = optimizer(constraints, objective, options, ...
    {YHCmpc{1}(:,1),YBmpc{1}(1),[pv_size,b_size]},{Eps{:}});
fprintf(' took %2f seconds \n', toc)

%%
if 1
% save the object
Fobj = saveobj(controller_obj);
Fstate = saveobj(controller_states);
Feps = saveobj(controller_eps);
save HF_optimizer_objects Fstate Fobj Feps 
end

