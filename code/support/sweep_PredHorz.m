%% Set-up
%select an arbitrary dsign and control parameters
d= [61.1, 409., 1.99, 0.457];
%select prediction horizon values
Np= [8,12,16,20,24];
%select wether data should be made or loaded
make_data = 0;
%number of evaluations for each value in the prediction horizon Np
Nreps = 5;
%to load run_OC add path 
addpath(genpath('../'))

%% Get data
if make_data
    objVal = zeros(length(Np), Nreps);
    evalTime = zeros(length(Np), Nreps);
    % Nreps evaluations of the design and control parameters under each Np
    for i = 1:length(Np)
       for j = 1:Nreps
          tic; disp(j)
          objVal(i,j) = run_OC(0,d,1,Np(i));
          evalTime(i,j) = toc;
       end
    end
    %average and STD of objective bvalues and solve times
    objVal_m = mean(objVal');
    evalTime_m = mean(evalTime');
    objVal_std = std(objVal');
    evalTime_std = std(evalTime');    
    
    save ../Output_Data/sweep_predHorz_data objVal_m evalTime_m objVal_std evalTime_std
else 
    
    load('../Output_Data/sweep_predHorz_data')
end

%% Plot

col1 = [0.4010 0.7450 0.9330];
col2 = [0.8500 0.3250 0.0980];

fig = figure; 
set(fig,'defaultAxesColorOrder',[col1;col2]);
yyaxis left
ylabel('simple regret')
hold on;
plot(Np,-(objVal_m(end)-objVal_m), 'color',col1,'LineWidth', 2);
errorbar(Np,-(objVal_m(end)-objVal_m),objVal_std,'color', col1,...
            'CapSize', 8, 'LineWidth', 1.25,'LineStyle','none');

yyaxis right
ylabel('normalized evaluation time')
plot(Np, evalTime_m/max(evalTime_m), 'color',col2,'LineWidth', 2);
errorbar(Np, evalTime_m/max(evalTime_m),evalTime_std/max(evalTime_m),'color', col2,...
            'CapSize', 8, 'LineWidth', 1.25,'LineStyle','none');
xlabel('prediction horizon \theta_p')
ylim([0.4 1.1])
xlim([7.5 24.5])
%set(gcf,'color','w');
box on;
set(gca,'FontSize',20)





