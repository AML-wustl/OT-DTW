function [init_id, init_center] = kMeans_Mallows_init(A, K)
% distance: Mallows distance (Wasserstein distance for continuous measure on R)
% K-means++ to calculate the initial seeding for K-means clustering
% Inputs:
% A: input time series samples, dim: num_sample*ts_len
% K: pre-defined number of clusters
% Outputs:
% init_id: best initialized cluster membership ids
% init_center: best initialized cluster centers

    [N,len] = size(A); % num-of-samples * sample_length
    centers = zeros(K,len); % save centers of each cluster
    D = zeros(N,K); % sample-to-center distance matrix
    s = RandStream('mlfg6331_64'); % Create the random number stream for reproducibility
    
    centers(1,:) = A(randsample(s,N,1),:);
    D(:,1) = OT_Mallows(centers(1,:), A);
    w = D(:,1) / sum(D(:,1));
    
    for k = 2:K
        centers(k,:) = A(randsample(s,N,1,true,w),:);
        D(:,k) = OT_Mallows(centers(k,:), A);
        [tmp_min, min_id] = min(D(:,1:k),[],2); 
        w = tmp_min / sum(tmp_min);
    end
    
    init_id = min_id;
    init_center= centers;
    
end