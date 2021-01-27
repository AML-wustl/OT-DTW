function [D_sqrd,T] = EROT_Wasserstein(a,centroids,M,ground_d,tol)
% distance: entropy_Wasserstein
% Inputs:
% a: sample sequence with dim: 1 * len_1.
% centroids: K cluster centroids with dim: K * len_2.
% ground_d: ground distance matrix
% eps: parameter in front of entropy
% tol: threshold for early stop
% Outputs:
% Dist:dim 1*K, return W_2^2
% T: transport matrix
    [~,len] = size(a);
    if nargin < 3
        % ground distance matrix
        [i,j] = meshgrid(1:len);
        ground_d = (i-j).^2; % (i-j).^2 or abs(i-j)
        ground_d = ground_d / max(ground_d(:));
        % barycenter parameters
        eps = 1e-2;
        tol = 1e-3;
        M = exp(-ground_d / eps);
    end
	[K,~] = size(centroids);
    D_sqrd = zeros(1,K);
    maxiter = 100; % maxiter
    T = zeros(len);
        
    for cluster_id = 1:K
        compt = 0;
        diff = inf; D_old = inf;
        b = centroids(cluster_id,:)';
        % intialize u and v
        v = ones(len,1);
        while (compt < maxiter) && (diff > tol)
            compt = compt + 1;
            T_prev = T;
            % update u and v
            u = a' ./ (M*v);
            v = b ./ (M'*u);
            % assemble T 
            T = diag(v) * M * diag(u); % transport matrix
%           % check diff between T and T_prev
%             diff = norm(T - T_prev,'fro');
            % Distance relative decrease
            D = sum(sum(ground_d.*T));
%             if norm(D./D_old-1,inf) < tol
%                 break;
%             end
%             D_old = D;
            % Marginal difference
            if norm(sum(abs(sum(T,1) - b')),inf) < tol
                break;
            end
        end
        D_sqrd(cluster_id) = D;
    end
    
end