%% solve outer
clear
clc
Nrepeat = 10;


%% run BO
%inputs
npar = 2                                   ; %number of bayesian parameters
deterministic = 1                          ; %bool value
mdl = @(vars) evalIDO(0,table2array(vars)) ; %model function name 
BO_algo = "bayesopt"                       ; %use bayesopt, turbo, or random
max_eval = 12                              ; %max objective evals
aqi_fun = 'expected-improvement'           ; %acquisition function
save_bool = 1                              ; %bool
plot_bool = 1                              ; %bool: plot convergance and range of all BO simulations

% define Optimization variables
names=["pv_size",  "b_size"];
bound_min = [0, eps ];
bound_max = [540, 1300]; 

opt_vars = [];
for ii = 1:npar
    opt_vars = [opt_vars optimizableVariable(names(ii), [bound_min(ii), bound_max(ii)],'Type','real')]; 
end

% BO loop
results = cell(Nrepeat,1);
objMin = zeros(max_eval,Nrepeat);
BO_time = tic;
for i = 1: Nrepeat
    rng(i*10)
    results{i} = bayesopt(mdl,opt_vars,...
        'Verbose',1,...
        'AcquisitionFunctionName',aqi_fun,... 
        'IsObjectiveDeterministic', deterministic,...  
        'MaxObjectiveEvaluations', max_eval,...
        'NumCoupledConstraints',0, ...
        'MaxTime', inf,...
        'NumSeedPoint',6,...
        'GPActiveSetSize', 300,...
        'PlotFcn', []);

    objMin(:,i) = results{i}.ObjectiveMinimumTrace;
end

toc(BO_time)
save ../Output_Data/IDO_EI_12

