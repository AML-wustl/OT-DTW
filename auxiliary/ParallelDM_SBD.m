function DM = ParallelDM_SBD(X)

[m, ~]=size(X);

DM = zeros(m,m);

delete(gcp('nocreate'));
parpool(6); % # of cores

for i = 1:m-1
    %disp(i);
    rowi = X(i,:);
    %for j=i+1:m
    parfor j = i+1:m
        rowj = X(j,:);
        DM(i,j) = 1 - max( NCCc(rowi,rowj) );
    end
end

delete(gcp('nocreate'));

for i = 1:m-1
    %disp(i);
    for j = i+1:m
        DM(j,i) = DM(i,j);
    end
end

end