%% Bayesian optimization of MPC %%
% runBO(model,deterministic, BO_algo, max_evals_per_BO, save_data, 
%       save_fx_name, save_x_name,num_par,plot_each_data, bound_min,bound_max)
%
%
% 
 function [] =runBO(model, deterministic,BO_algo,max_eval, save_bool, ...
                    converg_name, pname, npar,plot_bool,aqi_fun,...
                    bound_min,bound_max,opt_vars)

                
%User inputs
npar = 2                        ; %number of bayesian parameters
deterministic = 0                   ; %bool value
model = @(var) (evalIDO(0,var))      ; %model function name 
BO_algo = "bayesopt"                 ; %use bayesopt, turbo, or random
num_BO_simulations = 10             ; %number of BO simulations 
max_eval = 18              ; %max objective evals
aqi_fun = 'expected-improvement'    ;
save_bool = 0                       ; %bool
plot_bool = 1                   ; %bool: plot convergance and range of all BO simulations

%% Optimization variables

cont_names=['pv_size', 'b_size'];
bound_min_c = [0,eps];
bound_max_c = [540,1300];

% define Optimization variables
opt_vars = [];
for ii = 1:length(cont_names)
    opt_vars = [opt_vars optimizableVariable(cont_names(ii), [bound_min_c(ii), bound_max_c(ii)],'Type','real')]; 
end
%% Check User Inputs
% save objective evaluations to save_fx_name
% converg_name ="../data/"+BO_algo+"_"+num2str(deterministic)+'_'+...
        num2str(max_eval)+"_woUT-CSTR_fx.mat"   ; 
    
% save optimal input from BO to save_x_name    
pname ="../Output_Data/IDO"+BO_algo+"_"+num2str(deterministic)+'_'+...
             num2str(max_eval)+".mat"    ;

%%
% Perform bayesian optimization 
optimize = true;

BO_time = tic;
if optimize
    if isequal(BO_algo ,"bayesopt") || isequal(BO_algo,"random")
        
        if BO_algo == 'random'; nsp = max_eval;
        else; nsp = floor(max_eval/5); end

        det_cons = [];
        if deterministic; det_cons = [true]; 
        else; det_cons = [false]; end
        
        results = bayesopt(model,opt_vars,...
            'Verbose',0,...
            'AcquisitionFunctionName',aqi_fun,... %-plus',...
            'IsObjectiveDeterministic', deterministic,... % simulations with noise --> objective function is not deterministic
            'MaxObjectiveEvaluations', max_eval,...
            'NumCoupledConstraints',1, ...
            'AreCoupledConstraintsDeterministic', det_cons, ...
            'MaxTime', inf,...
            'NumSeedPoint',nsp,...
            'GPActiveSetSize', 300,...
            'PlotFcn', []);

             best_vars = results.bestPoint;

            p_star = table2array(best_vars);
            saveData(p_star, pname, 'ps')
            
            convergence = results.ObjectiveMinimumTrace';
            saveData(convergence, converg_name,'R' )
            
            if plot_bool;      plotBO('BayesOpt', 0, convergence);     end


    elseif BO_algo == "turbo"
        commandStr = 'python3.6 exec_fun_example.py ';
        systemCommand = [commandStr,num2str(bound_min'),' ',num2str(bound_max')]
        % the following line does not work, so matlab cant call teh python script 
         [status, results] = system(systemCommand);

        'Matlab is having an issue; please copy/paste the follwing command into a terminal directed to your turbo directory, where exec_fun.py is. Once TuRBO is done, continue the script. The data will be saved by python.  '
        systemCommand
        
        'leave a break-point here '

    end

    toc(BO_time)
    
end




