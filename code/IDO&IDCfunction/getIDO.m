function [mff, sff_hf] = getIDO()
%% Setup
addpath(genpath('../../'))
addpath(genpath('../../../mfBO'))

%% build problem
f_hf=@(d) evalIDO(0,d) ;
f_lf=@(d) evalIDO(1,d) ; 

%% mfBO setup
numFidels = 2;
bound_min = [0, eps];
bound_max = [540, 13*100]; 
funcHs = {f_lf; f_hf};   
hfMaxPt = [];
hfMaxVal = [];
bounds = [bound_min', bound_max'];
costs = [0.02; 1]; 

mff = mfFunction(funcHs, bounds, costs, [], hfMaxPt, hfMaxVal);
sff_hf = mfFunction({funcHs{numFidels}}, bounds, costs(numFidels), [], ...
                   hfMaxPt, hfMaxVal);

end