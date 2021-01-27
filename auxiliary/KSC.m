function [member_id, centers, clusters] = KSC(A, K)

m = size(A, 1);
member_id = ceil(K * rand(m, 1));
centers = zeros(K, size(A, 2));
D = zeros(m,K);

for iter = 1:100
    disp(iter);
    prev_mem = member_id;
    
    for k = 1:K
        centers(k,:) = ksc_centroid(member_id, A, k, centers(k,:));        
    end
    
    for i = 1:m
        %x = A(i,:);
        for k = 1:K
            %y = centers(k,:);
            dist = dhat_shift(centers(k,:),A(i,:));
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

function ksc = ksc_centroid(member_id, A, k, cur_center)
%Computes ksc centroid
a = [];
for i=1:length(member_id)
    if member_id(i) == k
        if sum(cur_center) == 0
            opt_a = A(i,:);
        else
            [tmp, tmps, opt_a] = dhat_shift(cur_center, A(i,:));
        end
        a = [a; opt_a];
    end
end
if size(a,1) == 0
    ksc = zeros(1, size(A,2)); 
    return;
end
b = a ./ repmat(sqrt(sum(a.^2,2)), [1 size(a,2)]);
M = b'*b - size(a, 1) * eye(size(a, 2));
[V, D] = eig(M);

ksc = V(:,end);

finddistance1 = sqrt(sum((a(1,:) - ksc').^2));
finddistance2 = sqrt(sum((a(1,:) - (-ksc')).^2));
   
if (finddistance1 >= finddistance2)
    ksc = -ksc;
end

%if sum(ksc) < 0
%    ksc = -ksc;
%end
end

function [dist, optshift, opty] = dhat_shift(x, y)
min_d = scale_d(x,y);
L = length(y);
range = -length(x):length(x);
for shift = range
    if shift < 0
        yshift = [y(-shift + 1:end) zeros(1, -shift)];
    else
        yshift = [zeros(1,shift) y(1:end-shift) ];
    end
    cur_d = scale_d(x,yshift);
    if cur_d <= min_d
        optshift = shift;
        opty = yshift;
        min_d = cur_d;
    end
end
dist = min_d;
end

function dist = scale_d(x,y)
alpha = x * y' / (y * y');
dist = norm(x - alpha * y) / norm(x);
end
