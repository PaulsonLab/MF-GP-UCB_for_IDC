function [mff, sff_hf] = getIDC(Var)
%% Setup
addpath(genpath('../../'))
addpath(genpath('../../../mfBO'))

%% build problem
f_hf=@(d) run_OC(0,d,Var) ;
f_lf=@(d) run_OC(1,d,Var) ; 

%% mfBO setup
numFidels = 2;
bound_min = [0, eps,0,0];
bound_max = [540, 13*100, 2, 2]; 
funcHs = {f_lf; f_hf};   
hfMaxPt = [];
hfMaxVal = [];
bounds = [bound_min', bound_max'];
costs = [0.04; 1]; 

mff = mfFunction(funcHs, bounds, costs, [], hfMaxPt, hfMaxVal);
sff_hf = mfFunction({funcHs{numFidels}}, bounds, costs(numFidels), [], ...
               hfMaxPt, hfMaxVal);


end