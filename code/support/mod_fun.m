%% function out = mod_fun(fun,x)    
%    all_outputs = fun(x);
%    out = all_outputs(1);
%
% simple function to supress unneeded outputs, primairly for blackbox
% optimization
%
function out = mod_fun(fun,x)
    
    all_outputs = fun(x);
    
    out = all_outputs(1);
    
end






