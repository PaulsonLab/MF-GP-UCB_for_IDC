%% function [] = Clustering(Nmax, n_st, n_end, save_bool)
% Nmax = number of clusters to be formed
% n_st = 356*24*start_yr, where start_year is in [1:4]
% n_end = 356*24*end_yr, where end_year is in [2:5]
function [Cluster_Lump,Stats,Centroids] = Clustering(Nmax, n_st, n_end )

% if isempty(nargin); save_bool = 0; 
% else; save_bool = nargin{1}; 
% end

load('../data/Uncertainty_Data.mat');
rng(100);

nyears = (n_end-(n_st))/(24);
%nyears = (n_end-(n_st-1))/(365*24);

% Set Uncertainty vectors into required size
AvgTemp = reshape(TempData(n_st:n_end-1,2),24,nyears);
DHI_Daily = reshape(DHI_Data(n_st:n_end-1,2),24,nyears);
DNI_Daily = reshape(DNI_Data(n_st:n_end-1,2),24,nyears);

% % Rescale vectors for better clustering
% AvgTemp = rescale(AvgTemp);
% DHI_Daily = rescale(DHI_Daily);
% DNI_Daily = rescale(DNI_Daily);

% Stack scaled uncertainties
Unc_St = [AvgTemp; DHI_Daily; DNI_Daily];

% Run k-means clustering
[idx,Centroids,sumd,D] = kmeans(Unc_St',Nmax);
SumD_Avg = mean(sumd);


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



k = 0;
order1 = ones(Nmax,1);
while Nmax~=length(unique(order1)) %make sure oders contains all clusters
    k=k+1;
    bin_sz = floor(length(idx)/Nmax/k); 
    for i = 1:Nmax*k
        order1(i) = mode(idx((i-1)*bin_sz+1:i*bin_sz));
    end
    if k>365*nyears
        disp('Excessive iterations when ordering ')
        break
    end
end
ordered_clusters = unique(order1,'stable');

prob = ones(Nmax,1);
for i = 1:Nmax % build cluster probabilities
    prob(i) = sum(idx == ordered_clusters(i))/length(idx);
end

figure(1); hold on;
scatter(1:nyears, idx, 'ob')
breaks = round(linspace(1,nyears,Nmax+1));
stairs([1;cumsum(prob*nyears)], [ordered_clusters;ordered_clusters(end)], '-r', 'LineWidth',2)
diff([1;cumsum(prob*nyears)]); %number of days in each cluster

for i = 1:Nmax
   CentroidsO(i,:)= Centroids(ordered_clusters(i),:); 
   Cluster_LumpO{i}=Cluster_Lump{ordered_clusters(i)};
end
Centroids = CentroidsO;
Cluster_Lump=Cluster_LumpO;

%Plot Centroids
figure(2);hold on;
for i = 1:Nmax
    cent_p = reshape(CentroidsO(i,:),24,3);
    for ii =1:3
        subplot(3,1,ii); hold on;
        plot(1:24, cent_p(:,ii)) 
    end
end



% %Unscale Centroids
% Tbnds   = [min(TempData(n_st:n_end,2)), max(TempData(n_st:n_end,2))];
% DHIbnds = [min(DHI_Data(n_st:n_end,2)), max(DHI_Data(n_st:n_end,2))];
% DNIbnds = [min(DNI_Data(n_st:n_end,2)), max(DNI_Data(n_st:n_end,2))];
% 
% Centroids(:,1:24)  = rescale(Centroids(:,1:24),Tbnds(1), Tbnds(2)); 
% Centroids(:,25:48) = rescale(Centroids(:,25:48),DHIbnds(1), DHIbnds(2));
% Centroids(:,49:72) = rescale(Centroids(:,49:72),DNIbnds(1), DNIbnds(2));

if 1
    save('../data/Clusterdata365_5','Cluster_Lump','Stats','Centroids');
end

