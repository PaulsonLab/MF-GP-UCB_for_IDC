% A script to test the MF-GP-UCB framework on IDO for the Solar Heating and Cooling System case study.

clear;
tic;
% randMethods = {'MF-GP-UCB', 'GP-UCB_HF', 'RAND'};
randMethods = {'MF-GP-UCB', 'GP-UCB_HF','RAND'};
numRandMethods = numel(randMethods);
Nrepeat = 10 ;

% Pick toy problem --------------------------------------------------------------
toyFuncHandle = @getIDO; toyFuncDesc = 'IDO';
% -------------------------------------------------------------------------------

addpath(genpath('../../'))
addpath(genpath('../../../mfBO'))

% Experiment parameters
[mff, sff_hf] = toyFuncHandle();
budget = 12;
params.initBudget = budget/3;
params.gpNoiseVars = [1E-4 1E-4];
numDims = mff.numDims;
params.budgetType = 'givenCost';
toyFuncDesc = sprintf('%s-%dD', toyFuncDesc, numDims);

% Set problem parameters
bounds = mff.bounds;
numFidels = mff.numFidels;
numDims = mff.numDims;
costs = mff.costs;

% To store the results for the random methods
history_mf = cell(Nrepeat,1);
sR_mf = cell(Nrepeat,1);
Cost_mf = cell(Nrepeat,1);
history_hf = cell(Nrepeat,1);
sR_hf = cell(Nrepeat,1);
Cost_hf = cell(Nrepeat,1);
history_rnd = cell(Nrepeat,1);
sR_rnd = cell(Nrepeat,1);
Cost_rnd = cell(Nrepeat,1);

fprintf(['(%s, max: %0.4f) \n', ...
  '=======================================================================\n'], ...
  toyFuncDesc, mff.hfMaxVal);

 
for methIter = 1:numRandMethods

  fprintf('\n');

  switch randMethods{methIter}
    
    case 'MF-GP-UCB'  
        
        for i = 1:Nrepeat
            rng(i*10);
            params.acqStrategy = 'MF-GP-UCB';
            mfFuncObj = mff;
            [maxVal, maxPt, queries, vals, history] = ...
             mfBO(mfFuncObj, bounds, budget, params);
            mfFuncObj = mff;
            history_mf{i} = history; 
            [sR, ~, ~, cC, ~] = getSimCumRegrets(mfFuncObj, history, params);
            sR_mf{i} = -sR(find(history_mf{i, 1}.evalFidels==2));
            Cost_mf{i} = cC(find(history_mf{i, 1}.evalFidels==2));
        end
    
    case 'GP-UCB_HF'
      
        for i = 1:Nrepeat 
          rng(i*10);
          params.acqStrategy = 'GP-UCB';
          mfFuncObj = sff_hf;
          [maxVal, maxPt, queries, vals, history] = ...
            mfBO(mfFuncObj, bounds, budget, params);
          history_hf{i} = history; 
          [sR, ~, ~, cC, ~] = getSimCumRegrets(mfFuncObj, history, params);
          sR_hf{i} = -sR;
          Cost_hf{i} = cC;
        end
        
    case 'RAND'
        
      for  i = 1:Nrepeat
          rng(i*10);
          [maxVal, maxPt, queries, vals, history] = ...
            randMaximise(sff_hf, bounds, budget, struct());
          mfFuncObj = sff_hf;
          history_rnd{i} = history;
          [sR, ~, ~, cC, ~] = getSimCumRegrets(mfFuncObj, history, params);
          sR_rnd{i} = -sR;
          Cost_rnd{i} = cC;
      end  
      
    otherwise
      error('Unknown Method');


  end % End Switch methods  


end % End for loop 

toc;
%% Plot the results out -------------------------------------------------------------------
save ../Output_Data/IDO_Full_12
% plot_sR('IDO_CaseStudy1')

