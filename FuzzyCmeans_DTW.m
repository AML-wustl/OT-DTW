function [init_id, ScoreMat, centers, obj_opt] = FuzzyCmeans_DTW(A, K, best_init_id)
% Input: A: dim N*len. Time series data.
% Input: K: pre-defined number of clusters.
% This function use fuzzy clustering to cluster input time series samples
% into K clusters, and record the scores belonging to each cluster.
% Outputs:
% init_id: current initialized cluster membership ids
% ScoreMatrix: dim N*K. The score (probability) of each sample
% centers: cluster centroids, dim: K * len
% obj_opt: sum of variance of all clusters 

    [N,len] = size(A); % num_of_samples * sample_length
    if nargin < 3
        member_id = ceil(K * rand(N, 1)); % save cluster id of each sample
    else
        member_id = best_init_id;
    end
    init_id = member_id;
    centers = zeros(K,len); % initialize centers of K clusters
    Dist = zeros(N,K); % ground distance matrix
    Max_iter = 100; % maximum iterations before stop
    window = len; % DTW window constraint
    obj = inf;
    m_exp = 2;
    tmp = rand(N,K);
    U = tmp ./ sum(tmp,2);
    eps = 1e-4;

    for iter = 1:Max_iter
        disp(iter);
        prev_U = U;
%         prev_obj = obj;

        for k = 1:K
            % Barycenter calculation
            centers(k,:) = DBA_fuzzy(A, U(:,k), m_exp);
        end

        for i = 1:N
            % calculate distance from each sample to each center
            [Dist(i,:),~] = cDTW(A(i,:), centers, window);
        end

        % Update U
        for i = 1:N
            for k = 1:K
                U(i,k) = 1 / sum((Dist(i,k)./Dist(i,:)).^(2/(m_exp-1)));
            end
        end
        
        obj = sum(sum(U.^m_exp .* (Dist.^2)));
        % Display intermediate outputs
        disp(obj);
        disp(norm(prev_U - U,'fro'));
        if norm(prev_U - U,'fro') <= N * K * eps
            break;
        end
    end
    
    % save the optimal objective
    obj_opt = obj;
    ScoreMat = U;
    
end