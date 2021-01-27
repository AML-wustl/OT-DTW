function [centroid] = Mallows_baryavg(X, weights)
% Computes centroid for Mallows barycenter
% INPUT:
% X: dim N*len, N histograms of size len
% OUTPUT:
% centroid: dim 1*len
    [n,len] = size(X);
    
    if n == 1
        centroid = X;
        return;
    end

    if nargin < 2
        weights = ones(n,1) / n; % uniform weights
    else
        weights = reshape(weights,n,1); % just to make sure we get a column weight vector.
    end
    
    % Construct unique sorted cumulative weights from all samples. 
    X_cum = cumsum(X,2);
    w_union = [0,X_cum(1,:)];
    for id = 2:n
        w_union = union(w_union, X_cum(id,:));
    end
    
    % Calculate the center and radius of each interval using inverse CDF.
    w_len = length(w_union);
    C = zeros(n, w_len-1);
    R = zeros(n, w_len-1);
    for id = 1:n
        for id_m = 1:(w_len-1)
            I_begin = invCDF(w_union(id_m),X_cum(id,:),len);
            I_end = invCDF(w_union(id_m+1),X_cum(id,:),len);
            C(id,id_m) = (I_begin + I_end) / 2;
            R(id,id_m) = (I_end - I_begin) / 2;
        end
    end
    
    % The C, R representation of the barycenter
    C_avg = mean(C,1); R_avg = mean(R,1);
    
    % build the centroid from C_avg and R_avg
    right_bound = C_avg + R_avg;
    centroid_cum = zeros(1,len+1);
    for len_id = 1:len
        if len_id <= right_bound(1)
            centroid_cum(len_id+1) = w_union(2) * len_id / right_bound(1);
        else
            union_id = binary_search(len_id, right_bound);
            if right_bound(union_id) ~= right_bound(union_id-1)
                centroid_cum(len_id+1) = w_union(union_id) + (w_union(union_id+1)-w_union(union_id)) * ...
                    (len_id-right_bound(union_id-1)) / (right_bound(union_id)-right_bound(union_id-1));
            else
                centroid_cum(len_id+1) = 1;
            end
        end
    end
    centroid = diff(centroid_cum);

end