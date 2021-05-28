function [] = process_history(H)
nreps = size(H,1);

for n =1:nreps
    bestval(n) = -H{n}.hfMaxVal;
    times(n)  = H{n}.evalTimes(end);
    
    values = -H{n}.evalVals;
    index = find(values==bestval(n));
    %bestpoints(n,:) = H{n}.evalPts(index,:);
    
    %evalIDO(0,bestpoints(n,:))
    %run_OC(0,bestpoints(n,:))
end

% mean solve time 
ci = std(bestval)/sqrt(nreps);

fprintf('Mean solve time for full Bayesian run : %2f \n', mean(times))
fprintf('Mean cost and CI for %2i Bayesian runs : %2f , +/- %2f  \n',nreps, mean(bestval), ci)


end




















