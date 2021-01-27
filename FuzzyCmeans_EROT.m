function [init_id, ScoreMat, centers, obj_opt] = FuzzyCmeans_EROT(A, K, best_init_id)
% input: A: dim N * len. Time series data.
% input: K: pre-defined number of clusters.
% This function use fuzzy clustering to cluster input time series samples
% into K clusters, and record the scores.
% return ScoreMatrix: dim N*K. The score (probability) of each sample
% belonging to each cluster.
    [N,len] = size(A); % num_of_samples * sample_length
    if nargin < 3
        member_id = ceil(K * rand(N, 1)); % save cluster id of each sample
    else
        member_id = best_init_id;
    end
    init_id = member_id;
    centers = zeros(K,len); % initialize centers of K clusters
    Dist = zeros(N,K); % ground distance matrix
    
    % Fuzzy coefficients
    m_exp = 2;
    tmp = rand(N,K);
    U = tmp ./ sum(tmp,2);
    eps = 1e-5;
    
    % OT coefficients
    [i,j] = meshgrid(1:len);
    ground_d = (i-j).^2; % (i-j).^2 or abs(i-j)
    ground_d = ground_d / max(ground_d(:)); % ground distance matrix
    % barycenter parameters
    Max_iter = 100;
    eps_OT = 1e-2;
    useGPU = false;
    tol = 1e-3;
    obj = inf;
    M = exp(-ground_d / eps_OT);
    
    % Main iteration loop
    for iter = 1:Max_iter
        disp(iter);
        prev_U = U;
%         prev_obj = obj;

        for k = 1:K
            centers(k,:) = EROT_baryavg(A(member_id==k,:),M,useGPU,tol);
        end

        for i = 1:N
            % calculate distance from each sample to each center
            [Dist(i,:),~] = EROT_Wasserstein(A(i,:), centers, M, ground_d, tol);
        end

        % Update U
        for i = 1:N
            for k = 1:K
                U(i,k) = 1 / sum((Dist(i,k)./Dist(i,:)).^(2/(m_exp-1)));
            end
        end
        
        obj = sum(sum(U.^m_exp .* (Dist.^2)));
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