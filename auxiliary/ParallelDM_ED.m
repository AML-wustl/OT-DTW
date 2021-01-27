function DM = ParallelDM_ED(X)

[m, ~] = size(X);

DM = zeros(m,m);

delete(gcp('nocreate'));
parpool(6); % # of cores

for i = 1:m-1
    %disp(i);
    rowi = X(i,:);
    %parfor j=i+1:m
    for j = i+1:m
        rowj = X(j,:);
        DM(i,j) = ED(rowi,rowj);
    end
end

 % Delete the current parallel pool if there is one, without creating one.
delete(gcp('nocreate'));

for i = 1:m-1
    %disp(i);
    for j = i+1:m
        DM(j,i) = DM(i,j);
    end
end

end