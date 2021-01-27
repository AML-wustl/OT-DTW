function [member_id, centers, clusters] = kMeans_ED(A, K)
% distance: Euclidean Distance
% time-series averaging: arithmetic mean

m = size(A, 1);
member_id = ceil(K * rand(m, 1));
centers = zeros(K, size(A, 2));
D = zeros(m,K);

for iter = 1:100
    disp(iter);
    prev_mem = member_id;
    
    for k = 1:K
        centers(k,:) = kMeans_ED_centroid(member_id, A, k);     
    end
    
    for i = 1:m
        for k = 1:K
            dist = ED(A(i,:),centers(k,:));
            D(i,k) = dist;
        end
    end
    
    [val, member_id] = min(D,[],2);
    if norm(prev_mem - member_id) == 0
        break;
    end
end

clusters = zeros(1,K);
for k = 1:K
    clusters(k) = length(find(member_id == k));
end

end

function centroid = kMeans_ED_centroid(member_id, A, k)
%Computes centroid
a = [];
for i = 1:length(member_id)
    if member_id(i) == k
        opt_a = A(i,:);
        a = [a; opt_a];
    end
end

if size(a,1) == 0
    centroid = zeros(1, size(A,2)); 
    return;
end

centroid = mean(a);

end
