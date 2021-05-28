%% solve outer
%% Setup
yhc_init =  [21;20;4]; %inside temp, inside wall temp, wall core temp
yb_init  = 50;

decomp_by_cluster = 1; % 1 for LF 0 for HF
sim_horz = 'y'; % y = year, m = month 
nCluster = 5;
%2,3,5,10,30 are the availible data sets for 1 month sim
%5, are the availible data sets for 1 year sim

addpath('../support')

%% build problem
if sim_horz == 'y'
    df_name = "../data/Clusterdata365_"+nCluster+".mat"; 
elseif sim_horz == 'm'
    df_name = "../data/Clusterdata"+nCluster+".mat"; 
end    
Q = 24;
[control_obj,control_states] = build_IDO_yalmip(decomp_by_cluster, df_name, nCluster, sim_horz);
% [control_obj_HF,control_states_HF] = build_IDO_yalmip(0, df_name, 365, sim_horz);
% [control_obj_LF,control_states_LF] = build_IDO_yalmip(1, df_name, 5, sim_horz);

%% run PSO
% f=@(d)control_obj{yhc_init, yb_init,d} ;
% DLB = [0,eps];
% DUB = [100, 13.5*(12*1000)*2];
% options = optimoptions('particleswarm','MaxIterations',50);
% tic
% fprintf('Solving IDO problem... ')
% [xmin,fmin,exitflag,output] = particleswarm(f,2,DLB,DUB,options);
% fprintf(' took %2f seconds \n', toc)
% fprintf('Minimum cost: %2f $ acheived with battery size: %2f kWh',fmin, xmin(2)/(12*1000))  
% fprintf(' and PV size: %2f \n', xmin(1))

%% run BO
%inputs
npar = 2                        ; %number of bayesian parameters
deterministic = 1                   ; %bool value
f=@(d) control_obj{yhc_init, yb_init,d} ;
mdl = @(vars) mod_fun(f,table2array(vars))      ; %model function name 
BO_algo = "bayesopt"                 ; %use baysopt, turbo, or random
max_eval = 100              ; %max objective evals
aqi_fun = 'expected-improvement'    ;
save_bool = 1                       ; %bool
plot_bool = 1                   ; %bool: plot convergance and range of all BO simulations

% Optimization variables
names=["pv_size",  "b_size"];
bound_min = [0, eps ];
bound_max = [540, 13*1000*12*10 ]; 

% define Optimization variables
opt_vars = [];
for ii = 1:npar
    opt_vars = [opt_vars optimizableVariable(names(ii), [bound_min(ii), bound_max(ii)],'Type','real')]; 
end

% save objective evaluations to save_fx_name
converg_name ="../data/IDO_"+sim_horz+"_"+BO_algo+"_nCluster"+num2str(nCluster)+"_J.mat"; 
% save optimal input from BO to save_x_name    
pname ="../data/IDO_"+sim_horz+"_"+BO_algo+"_nCluster"+num2str(nCluster)+"_x.mat"    ;


BO_time = tic;
% if optimize
%     if isequal(BO_algo ,"bayesopt") || isequal(BO_algo,"random")
%         
%         if BO_algo == 'random'; nsp = max_eval;
%         else; nsp = floor(max_eval/5); end
% 
%         det_cons = [];
%         if deterministic; det_cons = [true]; 
%         else; det_cons = [false]; end
        
        results = bayesopt(mdl,opt_vars,...
            'Verbose',1,...
            'AcquisitionFunctionName',aqi_fun,... %-plus',...
            'IsObjectiveDeterministic', deterministic,... % 
            'MaxObjectiveEvaluations', max_eval,...
            'NumCoupledConstraints',0, ...
            'MaxTime', inf,...
            'NumSeedPoint',10,...
            'GPActiveSetSize', 300,...
            'PlotFcn', []);

             best_vars = results.XAtMinObjective;

            p_star = table2array(best_vars);
            saveData(p_star, pname, 'ps')
            
            convergance = results.ObjectiveMinimumTrace';
            saveData(convergance, converg_name,'R' )
            
            if plot_bool;      plotBO('BayesOpt', 0, convergance);     end


    %end

    toc(BO_time)
    
%end


%% run best solution, collect data, and plot
fmin = convergance(end);
xmin = p_star;
fprintf('Minimum cost: %2f $ achieved with battery size: %2f kWh',fmin, xmin(2)/(12*1000))  
fprintf(' and PV size: %2f \n', xmin(1))

[sol] = control_states{yhc_init, yb_init,xmin};


X=zeros(4,nCluster*(Q+1));
X(1:3,:) = reshape(sol{1},3,nCluster*(Q+1));
X(4,:) = reshape(sol{2},1,nCluster*(Q+1));
U  = reshape(sol{3},3,Q,nCluster);
X = X([1,4],:);
if 1
    xhcLb = [19*ones(1,9),21*ones(1,10),19*ones(1,6)]; %
    xhcUb = [30*ones(1,9),26*ones(1,10),30*ones(1,6)]; %
    xbLb = ones(1,Q+1)*10;
    xbUb = ones(1,Q+1)*95;
    
    LB = repmat([xhcLb;xbLb],1,nCluster);
    UB = repmat([xhcUb;xbUb],1,nCluster);
    
    figure(1); hold on;
    for i = 1:2
        subplot(2,1,i); hold on;

        plot(1:(Q+1)*nCluster,X(i,:),'-or','linewidth',2)

        plot(1:(Q+1)*nCluster,LB(i,:),'k--','linewidth',2)
        plot(1:(Q+1)*nCluster,UB(i,:),'k--','linewidth',2)
        set(gcf,'color','w');
        set(gca,'FontSize',16)
        ylabel("X_"+i)
        xlabel('time [hrs]')
    end

%     UB = repmat([u_min{mod(1:Nsim,Q)+1}],1,Q);
%     LB = repmat([u_max{mod(1:Nsim,Q)+1}],1,Q);
    yl = [0 1500;0 1500;-1500 1500];
    figure(2); hold on;
    for i = 1:3
        subplot(3,1,i); hold on;
        plot(1:Q*nCluster,U(i,:),'-or','linewidth',2)
%         plot(1:Q,LB(i,:),'k--','linewidth',2)
%         plot(1:Q,UB(i,:),'k--','linewidth',2)
        set(gcf,'color','w');
        set(gca,'FontSize',16)
        ylim(yl(i,:))
        ylabel("u_"+i)
        xlabel('time [hrs]')
    end
end
