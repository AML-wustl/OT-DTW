function DM = ParallelDM_DTW(X,window)

[m, ~] = size(X);

DM = zeros(m,m);

delete(gcp('nocreate'));
parpool(6); % # of cores
%for i = 1:m
for i = 1:m-1
    disp(i);
    rowi = X(i,:);
    %for j=i+1:m
    parfor j = i+1:m
        rowj = X(j,:);
        DM(i,j) = cDTW(rowi,rowj,window);
    end
end

delete(gcp('nocreate'));
% Use the symmetry to fill the other half.
for i = 1:m-1
    %disp(i);
    for j = i+1:m
        DM(j,i) = DM(i,j);
    end
end

end