function [member_id, centers, clusters] = kMeans_LOROT(A, K)
% distance: LOROT_Wasserstein
% barycenter averaging: arithmetic mean

    [m,len] = size(A); % num of samples
    member_id = ceil(K * rand(m, 1)); % save cluster id of each sample
    centers = zeros(K,len); % save centers of each cluster
    D = zeros(m,K); 
    % ground distance matrix
    [i,j] = meshgrid(1:len);
    ground_d = abs(i-j); % (i-j).^2 or abs(i-j)
    % barycenter parameters
    Max_iter = 1000;
    lambda_1 = 1;
    lambda_2 = 24 / len;
    delta = 1;
    sigma = 2; % parameter for Gaussian distribution G
    useGPU = false;
    tol = 1e-6;

    for iter = 1:Max_iter
        disp(iter);
        prev_mem = member_id;

        for k = 1:K
%             centers(k,:) = meanavg(member_id, A, k);  
            tmp = LOROT_baryavg(A(member_id==k,:),ground_d,Max_iter,lambda_1,lambda_2,delta,sigma,useGPU,tol);
            centers(k,:) = tmp';
        end

        for i = 1:m
            [D(i,:),~] = LOROT_Wasserstein(A(i,:),centers,lambda_1,lambda_2,delta,sigma);
        end

        [~, member_id] = min(D,[],2);
        if norm(prev_mem - member_id) == 0
            break;
        end
    end
    % save how many samples in each cluster
    clusters = zeros(1,K); 
    for k = 1:K
        clusters(k) = length(find(member_id == k));
    end

end

function centroid = LOROT_baryavg(X,ground_d,Max_iter,lambda_1,lambda_2,delta,sigma,useGPU,tol,weights)
% Computes centroid for LOROT regularized barycenter
% INPUT:
% X: dim N * d , N histograms of size d
% ground_d: ground metric.
% Max_iter: maximum number of iterations
% lambda_1: 1st penalty parameter
% lambda_2: 2nd penalty parameter
% useGPU: 1, use; 0, no use.
% tol: stopping tolerance criterion
% weights: weights for weighted barycenter
%
% ADVICE: divide M by median(M) to have a natural scale
% for lambda
    n = size(X,1);

    if nargin < 10,
        weights = ones(n,1) / n;
    else
        weights = reshape(weights,n,1); % just to make sure we get a column vector.
    end
    MAX_COST = (24-1)^2;
    M = exp(-1/lambda_2 * (ground_d - lambda_1./(ground_d.^2/MAX_COST+delta))); % exp(-cost/gamma)
    G = 1/(sigma*sqrt(2*pi)) * exp(-ground_d/(2*sigma^2));
    K = G .* M; % the kernel projected on

    if useGPU    
        X = gpuArray(X);
        K = gpuArray(K);    
    end

    K(K<1e-300) = 1e-300;
    compt = 0;
    differ = inf;
    objectives = [];
    matrixVector = 0;
    MVP = [];

    % two first projections are simpler because u is necessarily a matrix of ones,
    % simplifying computations.
    T = K * (bsxfun(@rdivide,X',(sum(K))')); % first iteration, U = ones.
    u = bsxfun(@ldivide,T,exp(log(T)*weights));

    % loop below 
    while compt<Max_iter && differ>tol,
        % T: dim d*N
        T = u .* (K*(X'./(K'*u)));
        matrixVector = matrixVector + 2;
        compt = compt + 1;      
        u = bsxfun(@times,u,exp(log(T)*weights))./T;          
        matrixVector = matrixVector + 1;
        % log is also quite time-consuming. exp is only carried out on a vector
        MVP = [MVP,matrixVector];

        % what follows has only a marginal use and can be commented out.
        differ = sum(std(T,1,2)) % minimize the standard deviation of columns of T
        objectives = [objectives,differ];

        if mod(compt,5) == 1,
            differ = sum(std(T,1,2));    
            compt        
        end
    end
    centroid = mean(T,2); % dim: d*1

end