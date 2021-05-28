function [mff, sff_hf] = getIDC_3fidelity()
%% Setup
addpath(genpath('../../'))
addpath(genpath('../../../mfBO'))

%% build problem
f3=@(d) run_OC(0,d,1,24);
%f2=@(d) run_OC(0,d,1,24);     % clustered IDC (1)
f2=@(d) run_OC(1,d,1,12) ;      % no clusters IDO (2)
f1=@(d) run_OC(1,d,1,12) ;      % clustered IDO

%% mfBO setup
numFidels = 3;
bound_min = [0, eps,0,0];
bound_max = [540, 13*100, 2, 2]; 
funcHs = {f1;f2;f3};   
hfMaxPt = [];
hfMaxVal = [];
bounds = [bound_min', bound_max'];
%costs = [0.00998; .5794; 1]; %(1)
costs = [0.00998; .01637; 1]; %(2)

mff = mfFunction(funcHs, bounds, costs, [], hfMaxPt, hfMaxVal);
sff_hf = mfFunction({funcHs{numFidels}}, bounds, costs(numFidels), [], ...
               hfMaxPt, hfMaxVal);

end