function [Dist,D] = cDTW_kernel(a,centers,window)
% distance: constraint dynamic time warping (DTW)
% The warping matrix is determined from D.
    if nargin < 3
        window = 24;
    end
    [~,N] = size(a);
	[k,M] = size(centers);
    Dist = zeros(1,k);
    alpha = 2;

    for cluster_id = 1:k
        D = ones(N+1,M+1) * inf;
        D(1,1) = 0;
        for i = 2:N+1
            for j = max(2, i-window):min(M+1, i+window)
                p = a(i-1);
                q = centers(cluster_id,j-1);
%                 cost = (p - q)^2; % point-to-point local cost
                cost = (p - q)^2 + 2 * p * q * (abs(i - j)/24)^alpha; % point-to-point local cost
                D(i,j) = cost + min([D(i-1,j),D(i-1,j-1),D(i,j-1)]); % DP transition function
            end
        end
        Dist(cluster_id) = sqrt(D(N+1, M+1));
    end
    
end