% A script to test the MF-GP-UCB framework on synthetic functions.

clear;
tic;
% randMethods = {'MF-GP-UCB', 'GP-UCB_HF'};
randMethods = {'MF-GP-UCB'};
numRandMethods = numel(randMethods);
Nrepeat = 20;

FuncHandle = @getIDC_3fidelity; FuncDesc = 'IDC with three fidelities';
% -------------------------------------------------------------------------------

addpath(genpath('../../'))
addpath(genpath('../../../mfbo'))

% Experiment parameters
[mff, sff_hf] = FuncHandle();
budget = 30;
params.initBudget = 3;
numDims = mff.numDims;
params.budgetType = 'givenCost';
toyFuncDesc = sprintf('%s-%dD', FuncDesc, numDims);

% Set problem parameters
bounds = mff.bounds;
numFidels = mff.numFidels;
numDims = mff.numDims;
costs = mff.costs;

% To store the results for the random methods
simRegrets = cell(numRandMethods);
cumCosts = cell(numRandMethods);
history_mf = cell(Nrepeat,1);
maxObj_mf = cell(Nrepeat,1);
Cost_mf = cell(Nrepeat,1);
history_hf = cell(Nrepeat,1);
maxObj_hf = cell(Nrepeat,1);
Cost_hf = cell(Nrepeat,1);
fprintf(['(%s, max: %0.4f) \n', ...
  '=======================================================================\n'], ...
  FuncDesc, mff.hfMaxVal);

 
for methIter = 1:numRandMethods

  fprintf('\n');

  switch randMethods{methIter}
    
    case 'MF-GP-UCB'  
        
        for i = 1:Nrepeat
            params.acqStrategy = 'MF-GP-UCB';
            mfFuncObj = mff;
            [maxVal, maxPt, queries, vals, history] = ...
             mfBO( mfFuncObj, bounds, budget, params);
            mfFuncObj = mff;
            history_mf{i} = history; 
            [sR, ~, ~, cC, ~] = getSimCumRegrets(mfFuncObj, history, params);
            ev = history_mf{i}.evalVals(find(history_mf{i, 1}.evalFidels==3));
            maxval = zeros(length(ev),1);
            for j = 1:length(ev)
                maxval(j) = max(ev(1:j));   
            end    
            maxObj_mf{i} = -maxval;
            Cost_mf{i} = cC;
        end
    
       otherwise
        eror('Unknown Method');


  end % End Switch methods --------------------------------------------------------  


end % End for loop ----------------------------------------------------------------

toc;
%% Plot the results out -------------------------------------------------------------------

save ../Output_Data/IDC_CaseStudy3_2

%% post process
if 0
    load('../Output_Data/IDCplant_CaseStudy3_1')
    process_history(history_mf);
    C1 = Cost_mf;
    H1 = history_mf;
    obj1 = maxObj_mf;

    load('../Output_Data/IDCplant_CaseStudy3_2')
    process_history(history_mf);
    C2 = Cost_mf;
    H2 = history_mf;
    obj2 = maxObj_mf;

    figure(1); hold on
    for i = 1:10
        a = H1{i}.evalVals(find(H1{i}.evalFidels==3));
        b = H2{i}.evalVals(find(H2{i}.evalFidels==3));
        for j = 1:length(a)
            a(j) = max(a(1:j));
        end
        for j = 1:length(b)
            b(j) = max(b(1:j));
        end
        plot(C1{i}(find(H1{i}.evalFidels==3)), -a, '-r')
        plot(C2{i}(find(H2{i}.evalFidels==3)), -b, '-b')
    end

    figure(); hold on;
    plot_sR2(C1,H1,C2,H2,1);
    ylabel('SimpleRegret','FontSize',15);
    xlabel('Cost','FontSize',15);
    titleStr = sprintf('Simple Regret: Costs = %s', mat2str(mff.costs));
    title(titleStr);
    set(gcf,'color','w');
    set(gca,'yscale','log')
    box on;
    legend([p1 p3],{'clustered IDC','unclustered IDO'});
end

