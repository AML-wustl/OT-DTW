function [Dist,D] = cDTW(a,centroids,window)
% This function is to calculate cDTW distance between one sample and cluster centroids.
% Inputs:
% a: sample sequence with dim: 1 * len_1.
% centroids: K cluster centroids with dim: K * len_2.
% window: window constraint on time warping band. See Sakoe and Chiba 1978.
% Outputs:
% D: the warping matrix.
% Dist: the returned cDTW distance.
    
    [~,len_1] = size(a);
	[K,len_2] = size(centroids);
    Dist = zeros(1,K);

    for cluster_id = 1:K
        D = ones(len_1+1,len_2+1) * inf;
        D(1,1) = 0;
        for i = 2:len_1+1
            for j = max(2, i-window):min(len_2+1, i+window)
                p = a(i-1);
                q = centroids(cluster_id,j-1);
                cost = (p - q)^2; % point-to-point local cost
                D(i,j) = cost + min([D(i-1,j),D(i-1,j-1),D(i,j-1)]); % DP transition function
            end
        end
        Dist(cluster_id) = sqrt(D(len_1+1, len_2+1));
    end
    
end