function index = DB_index(member_id, centers, dat_input, clusters, M, ground_d, tol, dist_id)
% Iutputs:
% member_id: cluster membership ids, dim: N*1
% centers: cluster centers, dim: K*ts_len
% dat_input: input samples 
% clusters: saves the size of each cluster, dim: K*1
    [K,~] = size(centers);
    S = zeros(K,1);
    D = zeros(K,K);
    diff = zeros(K,K);
    if dist_id == 1 % OT distance
        for k = 1:K
            [Dist_sqrd,~] = EROT_Wasserstein(centers(k,:),dat_input(member_id == k,:),M,ground_d,tol);
            S(k) = sqrt(sum(Dist_sqrd) / clusters(k));
        end

        for i = 1:K
            for j = 1:i
                [D(i,j),~] = EROT_Wasserstein(centers(i,:),centers(j,:),M,ground_d,tol);
                D(j,i) = D(i,j);
                if j ~= i
                    diff(i,j) = (S(i) + S(j)) / sqrt(D(i,j));
                    diff(j,i) = diff(i,j);
                end
            end
        end
    elseif dist_id == 2 % DTW distance
        window = 24;
        for k = 1:K
            [Dist,~] = cDTW(centers(k,:),dat_input(member_id == k,:),window);
            S(k) = sum(Dist) / clusters(k);
        end

        for i = 1:K
            for j = 1:i
                [D(i,j),~] = cDTW(centers(i,:),centers(j,:),window);
                D(j,i) = D(i,j);
                if j ~= i
                    diff(i,j) = (S(i) + S(j)) / D(i,j);
                    diff(j,i) = diff(i,j);
                end
            end
        end
    else
        msg = 'Error occurred.';
        error(msg);
    end
    
    % calculate DB index
    R = max(diff,[],2);
    index = mean(R);
    
end