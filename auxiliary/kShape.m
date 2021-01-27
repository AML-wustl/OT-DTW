function [member_id, centers, clusters] = kShape(A, K)

m = size(A, 1);
member_id = ceil(K*rand(m, 1));
centers = zeros(K, size(A, 2));
D = zeros(m,K);

for iter = 1:100
    disp(iter);
    prev_mem = member_id;
    
    for k = 1:K
        centers(k,:) = kshape_centroid(member_id, A, k, centers(k,:));
    end
    
    for i = 1:m
        for k = 1:K
            dist = 1 - max( NCCc(A(i,:), centers(k,:) ) );
	    D(i,k) = dist;
        end
    end
    
    [val, member_id] = min(D,[],2);
    if norm(prev_mem-member_id) == 0
        break;
    end
end

clusters = zeros(1,K);
for k = 1:K
    clusters(k) = length(find(member_id == k));
end

end

function centroid = kshape_centroid(member_id, A, k, cur_center)
%Computes centroid
a = [];
for i=1:length(member_id)
    if member_id(i) == k
        if sum(cur_center) == 0
            opt_a = A(i,:);
        else
             [tmp, tmps, opt_a] = SBD(cur_center, A(i,:));
        end
        a = [a; opt_a];
    end
end

if size(a,1) == 0
    centroid = zeros(1, size(A,2)); 
    return;
end

[m, ncolumns]=size(a);
[Y, mean2, std2] = zscore(a,[],2);
S = Y' * Y;
P = (eye(ncolumns) - 1 / ncolumns * ones(ncolumns));
M = P*S*P;

[V, D] = qdwheig(M);

centroid = V(:,end);

finddistance1 = sqrt(sum((a(1,:) - centroid').^2));
finddistance2 = sqrt(sum((a(1,:) - (-centroid')).^2));

if (finddistance1 >= finddistance2)
    centroid = -centroid;
end

centroid = zscore(centroid);

end

