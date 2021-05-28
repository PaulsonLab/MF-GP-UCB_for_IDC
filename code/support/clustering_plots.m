clear
clc
addpath(genpath('../'))
load('Uncertainty_Data.mat');
Nmax = 5;
Nsample = 5;
rng(100);

% Set Uncertainty vectors into required size
AvgTemp = reshape(TempData(1:8760,2),24,365);
DHI_Daily = reshape(DHI_Data(1:8760,2),24,365);
DNI_Daily = reshape(DNI_Data(1:8760,2),24,365);

% Temp 365
figure; hold on;
for i = 1:365
    plot(1:24,AvgTemp(:,i),'-k','LineWidth',1.5)
end
xlabel('hour of the day');
ylabel('temperature [^oC]');
box on;
set(gcf,'Color','w');
set(gca, 'FontSize',20)

% DHI 365
figure; hold on;
for i = 1:365
    plot(1:24,DHI_Daily(:,i),'-k','LineWidth',0.1)
end
xlabel('hour of the day');
ylabel('DHI [W/m^2]');
box on;
alpha(0.25);
set(gcf,'Color','w');
set(gca, 'FontSize',20)

% DNI 365
figure; hold on;
for i = 1:365
    plot(1:24,DNI_Daily(:,i),'-k','LineWidth',0.1)
end
xlabel('hour of the day');
ylabel('DNI [W/m^2]');
box on;
set(gcf,'Color','w');
set(gca, 'FontSize',20)

%% K-means

% Stack scaled uncertainties
Unc_St = [AvgTemp; DHI_Daily; DNI_Daily];

% Run k-means clustering
SumD_Avg = zeros(Nmax,1);
for i = 1:Nmax
    [idx,C,sumd,D] = kmeans(Unc_St',i);
    SumD_Avg(i) = mean(sumd);
end

figure; hold on;
plot(1:Nmax,SumD_Avg/max(SumD_Avg),'-k','LineWidth',2);
xlabel('number of representative days');
ylabel('normalized sum of squared error');
box on;
set(gcf,'Color','w');
set(gca, 'FontSize',20)
% scatter(1:365,idx,'r');

% Collect all similar days together
Cluster_Lump = cell(Nmax,1);
for i = 1:Nmax
    Cluster_Lump{i} = Unc_St(:,idx(:) == i);
end

% Get mean and std of centroid distances
Centroid_Dist = [min(D')' idx];
Stats = zeros(Nmax,3);
for i = 1:Nmax
    ri = find(Centroid_Dist(:,2) == i);
    MCD = mean(Centroid_Dist(ri,1));
    SCD = std(Centroid_Dist(ri,1));
    Stats(i,:) = [i MCD SCD];
end
save('Clusterdata.mat','Cluster_Lump','Stats','C');

%% Clustered Plots
Rd = zeros(72,Nsample);
pd = zeros(1,Nsample);

for i = 1:Nsample
    [Rd(:,i), pd(i)] = InterclusterSample(i);
end


y = [16,C(1,8)-1.8; 17,C(2,8)-.4; 18,C(3,8)-1.1; 19,C(4,8)-1.3; 20,C(5,8)-1.3];
% Temp 5 Clusters
figure; hold on;
for i = 1:Nsample
    h = plot(C(i,1:24),'LineWidth',3);
    a = text(3.5,y(i,2), num2str(round(pd(i),3)));
    c = get(h,'Color');
    a.Color = c;
    a.FontSize = 14;
end
xlabel('hour of the day');
ylabel('temperature [^oC]');
box on;
% legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5');
set(gcf,'Color','w');
set(gca, 'FontSize',20)

% DHI 5 Clusters
figure; hold on;
for i = 1:Nsample
    plot(C(i,25:48),'LineWidth',3);
end
xlabel('hour of the day');
ylabel('DHI [W/m^2]');
box on;
% legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5')
set(gcf,'Color','w');
set(gca, 'FontSize',20)

% DNI 5 Clusters
figure; hold on;
for i = 1:Nsample
    plot(C(i,49:72),'LineWidth',3);
end
xlabel('hour of the day');
ylabel('DNI [W/m^2]');
box on;
% legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5')
set(gcf,'Color','w');
set(gca, 'FontSize',20)