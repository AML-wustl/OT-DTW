function Dist = ED(t,r)
% consider the case r could be multiples rows
	[r_r,~] = size(r);
    Dist = zeros(r_r,1);
    
    for sample_id = 1:r_r
        Dist(sample_id,:) = sqrt(sum((t - r(sample_id,:)).^2));
    end
end