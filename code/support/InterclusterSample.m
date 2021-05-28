%% [Rd,Pd] = InterclusterSample(CN)
% when CN in [1, length(Cluster_Lumps)]:
%       RD will be a random day from Cluster CN 
%       Pd will be the probability of a day from cluster CN occuring
% when CN = 0 :
%       RD will be a matrix, with each column i being a random day from
%       cluster i
%       Pd will be a vector with each element i being the prob. of a day
%       from cluster i occuring

function [Rd,Pd] = InterclusterSample(CN)
    
    load('../data/Clusterdata365_5.mat')
    
    if CN ~=0
        CL = Cluster_Lump{CN};
        sz = size(CL);
        for i = 1:sz(2)
            C_dist(i) = abs(norm(CL(:,i)) - norm(Centroids(CN,:)));
        end

        AugC = [C_dist; CL];
        AugC = sortrows(AugC')';
        Rd = AugC(2:end,ceil(abs(randn(1))*sz(2)/4));

        day_freq = size(Cluster_Lump{CN});
        Pd = day_freq(2)/365;
    else
        for ii = 1:length(Cluster_Lump)
            
            CN = ii;
            CL = Cluster_Lump{CN};
            sz = size(CL);
            C_dist = [];
            for i = 1:sz(2)
                C_dist(i) = abs(norm(CL(:,i)) - norm(Centroids(CN,:)));
            end
            
            AugC = [];
            AugC = [C_dist; CL];
            AugC = sortrows(AugC')';
            rd = AugC(2:end,ceil(abs(randn(1))*sz(2)/4));

            day_freq = size(Cluster_Lump{CN});
            pd = day_freq(2)/365;
            
            Rd(:,ii) = rd;
            Pd(:,ii) = pd;
            
        end
        
    end
    
end
