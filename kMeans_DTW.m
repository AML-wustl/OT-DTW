function [member_id, centers, obj_opt, clusters] = kMeans_DTW(A, K, best_init_id)
% distance: cDTW
% time-series averaging: DBA
% Inputs:
% A: input time series samples
% K: pre-defined number of clusters
% best_init_id: initialized cluster membership ids
% Outputs:
% init_id: current initialized cluster membership ids
% member_id: clustering output cluster membership ids
% centers: cluster centroids, dim: K * len
% obj_opt: sum of variance of all clusters 
% clusters: saves the size of each cluster

    [m,len] = size(A); % num_of_samples * sample_length
    if nargin < 3
        member_id = ceil(K * rand(m, 1)); % save cluster id of each sample
    else
        member_id = best_init_id;
    end
    centers = zeros(K,len); % initialize centers of K clusters
    Dist = zeros(m,K); % ground distance matrix
    Max_iter = 100; % maximum iterations before stop
    window = len; % DTW constraint
    obj = inf;

    for iter = 1:Max_iter
        disp(iter);
        prev_mem = member_id;
        prev_obj = obj;

        for k = 1:K
            % Barycenter calculation
%             centers(k,:) = meanavg(member_id, A, k); 
            centers(k,:) = DTW_baryavg(member_id, A, k);
        end

        for i = 1:m
            % calculate distance from each sample to each center
            [Dist(i,:),~] = cDTW(A(i,:),centers,window);
        end

        [dist_to_center, member_id] = min(Dist,[],2);
        obj = sumsqr(dist_to_center);
        disp(obj);
        if norm(prev_mem - member_id) == 0 || abs(obj-prev_obj) < 1e-3 * obj
            break;
        end
    end
    
    % save the optimal objective
    obj_opt = obj;
    % save how many samples in each cluster
    clusters = zeros(1,K);
    for k = 1:K
        clusters(k) = length(find(member_id == k));
    end

end

