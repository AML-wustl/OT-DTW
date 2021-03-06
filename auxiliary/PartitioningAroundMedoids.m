function [member_id, medoids, cost, clusters] = PartitioningAroundMedoids(medoids,D)
	n = size(D,1);
	K = size(medoids,1);
    maxit = 100;
	
	[costs,labels] = min(D(medoids,:));
	cost = sum(costs);
	
	last = 0;
	it = 0;
    while any(last ~= medoids & it < maxit)
		best_so_far_medoids = medoids;
		for i = 1:K
			medoids_aux = medoids;
			for j = 1:n
				medoids_aux(i) = j;
				[costs_aux,labels_aux] = min(D(medoids_aux,:));
				cost_aux = sum(costs_aux);
				if (cost_aux < cost)
					best_so_far_medoids = medoids_aux;
					cost = cost_aux;
					labels = labels_aux;
				end
			end
		end
		last = medoids;
		medoids = best_so_far_medoids;
		it = it + 1;
    end
    
    clusters = zeros(1,K);
    for k = 1:K
        clusters(k) = length(find(labels == k));
    end
    member_id = labels';
end