function [init_id, init_center] = kMeans_EROT_init(A, K)
% This code is the K-means++ to calculate the initial seeding for K-means 
% clustering under EROT_Wasserstein distance.
% Inputs:
% A: input time series samples, dim: num_sample*ts_len
% K: pre-defined number of clusters
% Outputs:
% init_id: best initialized cluster membership ids
% init_center: best initialized cluster centers

    [N,len] = size(A); % num-of-samples * sample_length
    centers = zeros(K,len); % save centers of each cluster
    D = zeros(N,K); % sample-to-center distance matrix
    % ground distance matrix
    [i,j] = meshgrid(1:len);
    ground_d = (i-j).^2; % (i-j).^2 or abs(i-j)
    ground_d = ground_d / max(ground_d(:));
    % barycenter parameters
    eps = 1e-2;
    tol = 1e-3;
    M = exp(-ground_d / eps);
    U = M .* ground_d;
    s = RandStream('mlfg6331_64'); % Create the random number stream for reproducibility
    
    centers(1,:) = A(randsample(s,N,1),:);
    [D(:,1),~] = EROT_Wasserstein(centers(1,:),A,M,ground_d,tol);
%     [D(:,1),~,~,~] = sinkhornTransportHotStartG(centers(1,:)',A',M,U);

    w = D(:,1) / sum(D(:,1));
    
    for k = 2:K
        centers(k,:) = A(randsample(s,N,1,true,w),:);
        [D(:,k),~] = EROT_Wasserstein(centers(k,:),A,M,ground_d,tol);
%         [D(:,k),~,~,~] = sinkhornTransportHotStartG(centers(k,:)',A',M,U);
        [tmp_min, min_id] = min(D(:,1:k),[],2); 
        w = tmp_min / sum(tmp_min);
    end
    
    init_id = min_id;
    init_center= centers;
    
end