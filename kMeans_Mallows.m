function [member_id, centers, obj_opt, clusters] = kMeans_Mallows(A, K, best_init_id)
% distance: Mallows distance (Wasserstein distance for continuous measure on R)
% barycenter averaging: closed-form expression
% Inputs:
% A: input time series samples, dim: num_sample*ts_len
% K: pre-defined number of clusters
% best_init_id: initialized cluster membership ids
% Outputs:
% init_id: current initialized cluster membership ids
% member_id: output cluster membership ids
% centers: cluster centroids, dim: K * len
% obj_opt: sum of variance of all clusters 
% clusters: saves the size of each cluster

    [N,len] = size(A); % num-of-samples * sample_length
    if nargin < 3
        member_id = ceil(K * rand(N, 1)); % save cluster id of each sample
    else
        member_id = best_init_id;
    end
    
    centers = zeros(K,len); % save centers of each cluster
    Dist_sqrd = zeros(N,K); % sample-to-center distance matrix

    % barycenter parameters
    Max_iter = 100;
    obj = inf;

    for iter = 1:Max_iter
        disp(iter);
        prev_mem = member_id;
        prev_obj = obj;

        for k = 1:K
            if size(find(member_id==k),1) == 0
                continue;
            else
                centers(k,:) = Mallows_baryavg(A(member_id==k,:));
            end
        end

        for i = 1:N
            % no transport matrix output for continuous measures
            Dist_sqrd(i,:) = OT_Mallows(A(i,:),centers);
        end

        [dist_to_center, member_id] = min(Dist_sqrd,[],2);
        obj = sum(dist_to_center);
        disp(obj);
        if norm(prev_mem - member_id) == 0 || abs(prev_obj - obj) < 1e-3 * obj
            break;
        end
    end
    % save the optimal objective
    obj_opt = sum(dist_to_center);
    
    % save how many samples in each cluster
    clusters = zeros(K,1);
    for k = 1:K
        clusters(k) = length(find(member_id == k));
    end

end
