function center = meanavg(member_id, A, k)
% Computes centroid
% no constraints on being original samples
    a = A(member_id==k,:);

    if size(a,1) == 0
        center = zeros(1, size(A,2)); 
        return;
    end

    center = mean(a);

end