% A script to test the MF-GP-UCB framework on synthetic functions.

clear;
start_time = tic;
% randMethods = {'MF-GP-UCB', 'GP-UCB_HF'};
randMethods = {'MF-GP-UCB', 'MF-GP-UCB'};
numRandMethods = numel(randMethods);
Nrepeat = 10;

FuncHandle = @getIDC; FuncDesc = 'IDC';
% -------------------------------------------------------------------------------

addpath(genpath('../../'))
addpath(genpath('../../../mfbo'))

 % To store the results for the random methods
simRegrets = cell(numRandMethods);
cumCosts = cell(numRandMethods);
history_mf = cell(Nrepeat,numRandMethods);
maxObj_mf = cell(Nrepeat,numRandMethods);
Cost_mf = cell(Nrepeat,numRandMethods);


for methIter = 1:2
    % Experiment parameters
    [mff, sff_hf] = FuncHandle(methIter);
    budget = 9;
    params.initBudget = budget/3;
    numDims = mff.numDims;
    params.budgetType = 'givenCost';
    toyFuncDesc = sprintf('%s-%dD', FuncDesc, numDims);

    % Set problem parameters
    bounds = mff.bounds;
    numFidels = mff.numFidels;
    numDims = mff.numDims;
    costs = mff.costs;

    fprintf(['(%s, max: %0.4f) \n', ...
      '=======================================================================\n'], ...
      FuncDesc, mff.hfMaxVal);


  fprintf('\n');

  switch randMethods{methIter}
    
    case 'MF-GP-UCB'  
        for i = 1:Nrepeat
            params.acqStrategy = 'MF-GP-UCB';
            mfFuncObj = mff;
            [maxVal, maxPt, queries, vals, history] = ...
             mfBO( mfFuncObj, bounds, budget, params);
            mfFuncObj = mff;
            history_mf{i,methIter} = history; 
            [sR, ~, ~, cC, ~] = getSimCumRegrets(mfFuncObj, history, params);
            ev = history_mf{i,methIter}.evalVals(find(history_mf{i, methIter}.evalFidels==2));
            maxval = zeros(length(ev),1);
            for j = 1:length(ev)
                maxval(j) = max(ev(1:j));   
            end   
            sR_mf{i,methIter} = - sR;
            maxObj_mf{i,methIter} = -maxval;
            Cost_mf{i,methIter} = cC;
        end
    
    otherwise
         error('Unknown Method');


  end % End Switch methods 

end % End for loop 

fprintf("Total time = %2i \n", toc(start_time));

%% save and plot
save ../Output_Data/IDC_CaseStudy2(2)

for n = 1:Nrepeat
    for i=1:2
        bestval(n,i) = history_mf{n,i}.hfMaxVal;

        values = history_mf{n,i}.evalVals;
        index = find(values==bestval(n,i));
        bestpoints = history_mf{n,i}.evalPts(index,:);
        bat(n,i) = bestpoints(2);
        pv(n,i) = bestpoints(1);
        bkl(n,i) = bestpoints(3);
        bku(n,i) = bestpoints(4);
    end
end

% Compare Designs and cost
figure(); hold on;
scatter(bat(:,1), bestval(:,1),250,'b','filled', 'linewidth',2)
scatter(bat(:,2), bestval(:,2),250,'r','filled', 'linewidth',2)
legend('Low Varience', 'High Varience')
ylabel('Cost')
xlabel('Battery size [kWh]')
box on;
set(gcf,'color','w');
set(gca,'FontSize',16)

figure(); hold on;
scatter(pv(:,1), bestval(:,1),250,'b','filled', 'linewidth',2)
scatter(pv(:,2), bestval(:,2),250,'r','filled', 'linewidth',2)
legend('Low Varience', 'High Varience')
ylabel('Cost')
xlabel('PV array size [M^2]')
box on;
set(gcf,'color','w');
set(gca,'FontSize',16)
        
figure(); hold on;
scatter(bat(:,1), pv(:,1),250,'b','filled', 'linewidth',2)
scatter(bat(:,2), pv(:,2),250,'r','filled', 'linewidth',2)
legend('Low Varience', 'High Varience')
xlabel('Battery size [kWh]')
ylabel('PV array size [M^2]')
box on;
set(gcf,'color','w');
set(gca,'FontSize',16)
        
% Compare Designs and cost probabilities
figure(); hold on;
[f,xi] = ksdensity(pv(:,1),'support', [0,540],'BoundaryCorrection', 'Reflection');
plot(xi,f,'b', 'linewidth',2)
[f,xi] = ksdensity(pv(:,2),'support', [0,540],'BoundaryCorrection', 'Reflection');
plot(xi,f,'r', 'linewidth',2)
xlabel('PV array size')
ylabel('Probability')
legend('Low Varience', 'High Varience')
box on;
set(gcf,'color','w');
set(gca,'FontSize',16)

figure(); hold on;
[f,xi] = ksdensity(bat(:,1),'support', [0,1300],'BoundaryCorrection', 'Reflection');
plot(xi,f, 'b', 'linewidth',2)
[f,xi] = ksdensity(bat(:,2),'support', [0,1300],'BoundaryCorrection', 'Reflection');
plot(xi,f, 'r', 'linewidth',2)
xlabel('Battery size')
ylabel('Probability')
legend('Low Varience', 'High Varience')
box on;
set(gcf,'color','w');
set(gca,'FontSize',16)

figure(); hold on;
[f,xi] = ksdensity(bestval(:,1));
plot(xi,f,'b', 'linewidth',2)
[f,xi] = ksdensity(bestval(:,2));
plot(xi,f,'r', 'linewidth',2)
xlabel('Cost')
ylabel('Probability')
legend('Low Varience', 'High Varience')
box on;
set(gcf,'color','w');
set(gca,'FontSize',16)