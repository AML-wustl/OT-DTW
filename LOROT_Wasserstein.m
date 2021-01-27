function [Dist,T] = LOROT_Wasserstein(a,centers,lambda_1,lambda_2,delta,sigma)
% distance: LOROT_Wasserstein
    
    [~,N] = size(a);
	[k,~] = size(centers);
    Dist = zeros(1,k);
    
    maxiter = 100; % maxiter

    for cluster_id = 1:k
        [i,j] = meshgrid(1:N);
        ground_d = abs(i-j); % ground distance
        MAX_COST = (24-1)^2;
        M = exp(-1/lambda_2 * (ground_d - lambda_1./(ground_d.^2/MAX_COST+delta))); % exp(-cost/gamma)
        G = 1/(sigma*sqrt(2*pi)) * exp(-ground_d/(2*sigma^2));
        % intialize u and v
        v = exp(-0.5/lambda_2) * ones(N,1);
        K = G .* M; % the kernel projected on
        for k = 1:maxiter
            % update u and v
            u = a'./(K*v);
            v = centers(cluster_id,:)'./(K'*u);
            % assemble pi
            T = diag(u) * K * diag(v);
%             % display pi (with marginals on top and to the left)
%             imagesc([p'/max(p) 0;pi/max(pi(:)) q/max(q)])
%             colormap(1-map)
%             drawnow
%             pause(0.1);
        end
        Dist(cluster_id) = sum(sum(ground_d.*T));
    end
    
end