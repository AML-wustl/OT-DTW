function [centroid] = EROT_baryavg(X,M,useGPU,tol,weights)
% Computes centroid for entropy regularized barycenter
% Modified From Peyre's github code
% INPUT:
% X: dim N * d , N histograms of size d
% ground_d: ground metric.
% Max_iter: maximum number of iterations
% eps: penalty parameter
% useGPU: 1, use; 0, no use.
% tol: stopping tolerance criterion
% weights: weights for weighted barycenter
%
% ADVICE: divide M by median(M) to have a natural scale
% for lambda
    n = size(X,1);

    if nargin < 6
        weights = ones(n,1) / n;
    else
        weights = reshape(weights,n,1); % just to make sure we get a column vector.
    end
    
    if useGPU    
        X = gpuArray(X);
        M = gpuArray(M);    
    end

    M(M<1e-300) = 1e-300;
    compt = 0;
    diff = inf;
    Max_iter = 1000;

    % two first projections are simpler because u is necessarily a matrix of ones,
    % simplifying computations.
    T = M * (bsxfun(@rdivide,X',(sum(M))')); % first iteration, U = ones.
    u = bsxfun(@ldivide,T,exp(log(T)*weights));

    % loop below 
    while (compt < Max_iter) && (diff > tol)
        % T: dim d*N.
        T = u .* (M*(X'./(M'*u)));
        u = bsxfun(@times,u,exp(log(T)*weights))./T;          
        % log is also quite time-consuming. exp is only carried out on a vector
        % what follows has only a marginal use and can be commented out.
        diff = sum(std(T,1,2)); % minimize the variance of columns of T
%         disp(diff);
        compt = compt + 1;
    end
    tmp = mean(T,2); % dim: d*1
    centroid = tmp';

end