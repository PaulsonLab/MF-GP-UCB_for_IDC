function [Fctrl,Fstate] = build_OC(n_pred)
import casadi.*
addpath(genpath('../'))

persistent idx_sim;
if isempty(idx_sim)
    idx_sim = 1;
end

%% Setup
% Periodic time horizon
Q = 24;

% Sampling time
Ts = 1;

% system dimensions
nx = 3;
nu = 3;
ny = 1;

% conversion ah/kwh
ah_kwh = (1000/12);

%% Formulate and simulate economic MPC problem
% Problem horizon
N = n_pred;

% Define optimal control problem (with periodic constraints)
x = sdpvar(repmat(3,1,N+1),repmat(1,1,N+1));
xb = sdpvar(repmat(1,1,N+1),repmat(1,1,N+1));
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
Eps = sdpvar(repmat(4,1,N),repmat(1,1,N));

xlb = sdpvar(repmat(1,1,N+1),repmat(1,1,N+1));
xub = sdpvar(repmat(1,1,N+1),repmat(1,1,N+1));
xblb = 10;
xbub = 95;

dni = sdpvar(repmat(1,1,N),repmat(1,1,N));
dhi = sdpvar(repmat(1,1,N),repmat(1,1,N));
human = sdpvar(repmat(1,1,N),repmat(1,1,N));
Tod = sdpvar(repmat(1,1,N),repmat(1,1,N));
cost = sdpvar(repmat(1,1,N),repmat(1,1,N));
Wxhc = sdpvar(24,1);

pv_size = sdpvar(1);
b_size = sdpvar(1);

b_coef = 100/(b_size*ah_kwh);

% Input constraints
ulb = [0;0;-1500];
uub = [1500;1500;1500] ;

constraints = [];
constraints = [constraints, 0<=[Eps{:}]];

objective   = 0; 

C1 = 9.356*10^5; C2 = 2.970*10^6; C3 = 6.695*10^5;
K1 = 16.48; K2 = 108.5; K3 = 5; K4 = 30.5; K5 = 23.04;
eta_c = 2;
eta_h = 4;

for k = 1:n_pred
    wk = [Tod{k};dni{k};human{k};dhi{k}];
    
    constraints = [constraints, ulb <= u{k}(:) <=uub];
    
    %Temperature continuity constraints
    constraints = [constraints,x{k+1}(1) == ...    
            3600*(1/C1)*((K1+K2)*(x{k}(2)-x{k}(1)) + K3*(wk(1)-x{k}(1))+...
            K5*(x{k}(3)-x{k}(1))+wk(2)+u{k}(1)*eta_h-u{k}(2)*eta_c + wk(3))...
            + x{k}(1)];
    constraints = [constraints,x{k+1}(2) == ...
            3600*(1/C2)*((K1+K2)*(x{k}(1)-x{k}(2)) + wk(2)) + x{k}(2)]; 
    constraints = [constraints,x{k+1}(3) == ...
            3600*(1/C3)*(K5*(x{k}(1)-x{k}(3))+K4*(wk(1)-x{k}(3))) + x{k}(3)];     
    %Temp state constraints  
    constraints = [constraints, (xlb{k+1}-Eps{k}(1) <= ...
           x{k+1}(1) <= xub{k+1}+Eps{k}(2) )];

    %SOC continuity constraints
    constraints = [constraints,xb{k+1} == ... %continuity
           xb{k}+b_coef*(-(u{k}(1)+u{k}(2)) +...
           u{k}(3) + wk(4)*pv_size*0.18)*ah_kwh];
    %SOC state constraints   
    constraints = [constraints, (xblb <= xb{k+1}...
            <= xbub+Eps{k}(4)  )]; % infesible otherwise
            
    objective = objective+(cost{k}*u{k}(3))+ [Wxhc(k),Wxhc(k)]*Eps{k}(1:2)+ Eps{k}(3:4)'*diag([2E6,2E6])*Eps{k}(3:4);

end

options = sdpsettings('solver','gurobi','verbose',0,'debug',1,'warning',1);
Fcontroller = optimizer(constraints, objective, options, {[[x{1}];[xb{1}]], [xlb{:}],...
    [xub{:}],Wxhc(1:n_pred), [dni{:}],[dhi{:}],[human{:}], ...
    [Tod{:}],[cost{:}],[pv_size, b_size ]}, [u{:}]);
Fpredictor = optimizer(constraints, objective, options, {[[x{1}];[xb{1}]], [xlb{:}],...
    [xub{:}],Wxhc, [dni{:}],[dhi{:}],[human{:}], ...
    [Tod{:}],[cost{:}],[pv_size, b_size ]}, [x{:};xb{:}]);



% save the object
Fctrl = saveobj(Fcontroller);
Fstate = saveobj(Fpredictor);
%save optimal_control_object_shortpredictions Fctrl Fstate 
Fctrl = (Fcontroller);
Fstate = (Fpredictor);


end