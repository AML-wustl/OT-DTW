function index = CH_index(member_id, center_spread, obj_opt, clusters)
% Iutputs:
% member_id: cluster membership ids, dim: N*1
% center_spread: W_2^2(centers(k,:), center_all), dim: K*1
% obj_opt: sum of variance of all clusters 
% clusters: saves the size of each cluster, dim: K*1
    N = length(member_id);
    [K,~] = size(center_spread);
    % BI: between cluster variation
    BI = 0;
    for k = 1:K
        BI = BI + clusters(k) * center_spread(k);
    end
    % WI: within cluster variation
    WI = obj_opt;
    if K ~= 1
        index = (N-K) / (K-1) * BI / WI;
    else
        error('only one cluster, CH index not valid');
    end
    
end