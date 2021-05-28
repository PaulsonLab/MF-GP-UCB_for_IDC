%%  cost = evalIDO( decomp, d)
% Function Caller
% Used to evaluate either high or low fidelity IDO problem
% decomp = 0,1 - [high fidelity, low fidelity]
% d = [pv array size [m^2], Battery size [kwh]]
% cost = annual cost of PV and battery installation
%

function cost = evalIDO( decomp, d)
addpath(genpath('../../'))
load_fstate = 1;
X0 = [21; 20; 4; 50];
if decomp % build low fidelity model (5 clusters)
    df = "Clusterdata365_5.mat";
    if load_fstate
        load('../optimizerObjects/LF_optimizer_objects.mat') %
    else
       Fstate = build_IDO_yalmip(1, df, 5, 'y'); 
    end
else %build High fidelity model
    df = "Clusterdata365_365.mat";
    if load_fstate
        load('../optimizerObjects/HF_optimizer_objects.mat')
    else
        Fstate = build_IDO_yalmip(0, df, 365, 'y');
    end
end

%determine cost by simulating physical constraint satisfaction
sol = Fstate(X0(1:3),X0(4),d(1:2));
%return cost
cost = ControlSimulator(sol,d(1:2));
cost = -cost;

end
