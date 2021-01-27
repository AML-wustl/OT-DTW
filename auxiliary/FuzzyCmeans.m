function [Member_Deg, centers] = FuzzyCmeans(X, K)
% input X: raw data
% input K: number of clusters
    [num_sample, length] = size(X);
    U_next = zeros(num_sample, K); % save degree of membership in each row
    % initialize U_next using random value generator
    for idx = 1:num_sample
        U_next(idx,:) = rand(1,K);
        U_next(idx,:) = U_next(idx,:) / sum(U_next(idx,:));
    end
    U =  inf * ones(num_sample, K);
    centers = zeros(K,length);
    m = 1.5;
    eps = 10^-5;
    while norm(U_next - U,'fro') > eps
        disp(U_next);
        U = U_next;
        % At k step, calculate the center vectors
        for id = 1:K
            centers(id,:) = (U(:,id).^m)' * X ./ sum(U(:,id).^m);
        end
        % Update U_next
        for i = 1:num_sample
            denom_ik = sum(cDTW(X(i,:),centers,length).^(2/(1-m)));
%             denom_ik = sum(ED(X(i,:),centers).^(2/(1-m)));
            for j = 1:K
                U_next(i,j) = cDTW(X(i,:),centers(j,:),length)^(2/(1-m)) / denom_ik;
%                 U_next(i,j) = ED(X(i,:),centers(j,:))^(2/(1-m)) / denom_ik;
            end
        end
    end
    Member_Deg = U_next;
end