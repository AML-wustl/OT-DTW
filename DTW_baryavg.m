function ksc = DTW_baryavg(member_id, A, k)
% This function is a connection between kMeans_DTW and DBA
% Outputs:
% ksc: returned DTW barycenter
    a = A(member_id==k,:);

    if size(a,1) == 0
        ksc = zeros(1, size(A,2)); 
        return;
    end

    ksc = DBA(a);

end


