clear
% load variables from corresponding data files
load('../Output_Data/IDO_CaseStudy1.mat','Nrepeat','history_hf');
load('../Output_Data/IDO_CaseStudy_EI2.mat','results');
load('../Output_Data/IDO_Rand_CaseStudy.mat','history_rnd');
iter = 18;
psize = zeros(3*iter,Nrepeat);
bsize = zeros(3*iter,Nrepeat);
cost = zeros(3*iter,Nrepeat);

for i = 1:Nrepeat
    psize(1:iter,i) = history_hf{i}.evalPts(:,1);
    psize(iter+1:iter*2,i) = history_rnd{i}.evalPts(:,1);
    psize((2*iter)+1:3*iter,i) = table2array(results{i}.XTrace(:,1));
    bsize(1:iter,i) = history_hf{i}.evalPts(:,2);
    bsize(iter+1:iter*2,i) = history_rnd{i}.evalPts(:,2);
    bsize((2*iter)+1:3*iter,i) = table2array(results{i}.XTrace(:,2));
    cost(1:iter,i) = history_hf{i}.evalVals;
    cost(iter+1:iter*2,i) = history_rnd{i}.evalVals;
    cost((2*iter)+1:3*iter,i) = results{i}.ObjectiveTrace;
end
psize = reshape(psize,540,1);
bsize = reshape(bsize,540,1);
cost = reshape(cost,540,1);
pv = 540.*rand(200000,1);
b = 1300.*rand(200000,1);
newpt = [pv b];
F = scatteredInterpolant([psize, bsize], cost,'linear','linear');
z = F(newpt);
figure(1)
scatter(newpt(:,1),newpt(:,2),2,z);
hold on
scatter(psize,bsize,'r','filled')
colorbar
title('Topology for original cost structure')
xlabel('PV Size');
ylabel('Battery Size');

%% Rand Scatterplot
clear
load('../Output_Data/IDO_Rand_TopologyStudy.mat','history_rnd');
Psize = zeros(300,1);
Bsize = zeros(300,1);
Cost = zeros(300,1);
Psize = history_rnd{1}.evalPts(:,1);
Bsize = history_rnd{1}.evalPts(:,2);
Cost = history_rnd{1}.evalVals;
pv = 540.*rand(200000,1);
b = 1300.*rand(200000,1);
newpt = [pv b];
F = scatteredInterpolant([Psize, Bsize], Cost,'linear','linear');
z = F(newpt);
figure(2)
scatter(newpt(:,1),newpt(:,2),2,z);
hold on
scatter(Psize,Bsize,'r','filled')
colorbar
title('Topology for updated cost structure')
xlabel('PV Size');
ylabel('Battery Size');

%% Rand Scatterplot 2
clear
load('../Output_Data/IDO_RAND_UpdatedTopology2.mat','history_rnd');
Psize = zeros(399,1);
Bsize = zeros(399,1);
Cost = zeros(399,1);
Psize = history_rnd{1}.evalPts(:,1);
Bsize = history_rnd{1}.evalPts(:,2);
Cost = history_rnd{1}.evalVals;
pv = 540.*rand(200000,1);
b = 1300.*rand(200000,1);
newpt = [pv b];
F = scatteredInterpolant([Psize, Bsize], Cost,'linear','linear');
z = F(newpt);
figure(3)
scatter(newpt(:,1),newpt(:,2),2,z);
hold on
scatter(Psize,Bsize,'r','filled')
colorbar
title('Topology for final cost structure')
xlabel('PV Size');
ylabel('Battery Size');


%% Rand Scatterplot 3
clear
load('../Output_Data/IDO_RAND_LFUpdatedTopology2.mat','history_rnd');
Psize = zeros(501,1);
Bsize = zeros(501,1);
Cost = zeros(501,1);
Psize = history_rnd{1}.evalPts(:,1);
Bsize = history_rnd{1}.evalPts(:,2);
Cost = history_rnd{1}.evalVals;
pv = 540.*rand(200000,1);
b = 1300.*rand(200000,1);
newpt = [pv b];
F = scatteredInterpolant([Psize, Bsize], Cost,'linear','linear');
z = F(newpt);
figure(4)
scatter(newpt(:,1),newpt(:,2),2,z);
hold on
scatter(Psize,Bsize,'r','filled')
colorbar
title('Topology for final cost structure LF')
xlabel('PV Size');
ylabel('Battery Size');


%% Rand Scatterplot 4
clear
load('../Output_Data/IDO_RAND_LFTopology.mat','history_rnd');
Psize = zeros(501,1);
Bsize = zeros(501,1);
Cost = zeros(501,1);
Psize = history_rnd{1}.evalPts(:,1);
Bsize = history_rnd{1}.evalPts(:,2);
Cost = history_rnd{1}.evalVals;
pv = 540.*rand(200000,1);
b = 1300.*rand(200000,1);
newpt = [pv b];
F = scatteredInterpolant([Psize, Bsize], Cost,'linear','linear');
z = F(newpt);
figure(5)
scatter(newpt(:,1),newpt(:,2),2,z);
hold on
scatter(Psize,Bsize,'r','filled')
colorbar
title('Topology for original cost structure LF')
xlabel('PV Size');
ylabel('Battery Size');
