%*******************************************************************************
 % This DBA copy is an extension of DBA.m to fuzzy weighted version
 %*****************************************************************************

function average = DBA_fuzzy(sequences,u,m)
	[~,index] = max(u);
    % choose the maximum index u_ik as the initial barycenter
	average = sequences(index,:);
	for i = 1:15
		average = DBA_one_iteration(average,sequences,u,m);
	end
end


function average = DBA_one_iteration(averageS,sequences,u,m)

	tupleAssociation = cell (1, size(averageS,2));
	for t = 1:size(averageS,2)
		tupleAssociation{t} = [];
	end

	costMatrix = [];
	pathMatrix = [];

	for k = 1:size(sequences,1)
	    sequence = sequences(k,:);
	    costMatrix(1,1) = distanceTo(averageS(1),sequence(1));
	    pathMatrix(1,1) = -1;
	    for i = 2:size(averageS,2)
            costMatrix(i,1) = costMatrix(i-1,1) + distanceTo(averageS(i),sequence(1));
            pathMatrix(i,1) = 2;
	    end
	    
	    for j = 2:size(sequence,2)
            costMatrix(1,j) = costMatrix(1,j-1) + distanceTo(sequence(j),averageS(1));
            pathMatrix(1,j) = 1;
	    end
	    
	    for i = 2:size(averageS,2)
            for j = 2:size(sequence,2)
                indiceRes = ArgMin3(costMatrix(i-1,j-1),costMatrix(i,j-1),costMatrix(i-1,j));
                pathMatrix(i,j) = indiceRes;

                if indiceRes == 0
                res = costMatrix(i-1,j-1);
                elseif indiceRes == 1
                res = costMatrix(i,j-1);
                elseif indiceRes == 2
                res = costMatrix(i-1,j);
                end

                costMatrix(i,j) = res + distanceTo(averageS(i),sequence(j));
            end
	    end
	    
	    i = size(averageS,2);
	    j = size(sequence,2);
	    
	    while(true)
            tupleAssociation{i}(end+1) = u(k)^m * sequence(j);
            if pathMatrix(i,j) == 0
                i=i-1;
                j=j-1;
            elseif pathMatrix(i,j) == 1
                j=j-1;
            elseif pathMatrix(i,j) == 2
                i=i-1;          
            else
                break;
            end
        end

	end

	for t = 1:size(averageS,2)
	   averageS(t) = sum(tupleAssociation{t}) / sum(u.^m);
	end
	   
	average = averageS;

end

function value = ArgMin3(a,b,c)

	if (a<b)
	    if (a<c)
		value=0;
		return;
	    else
		value=2;
		return;
	    end
	else
	    if (b<c)
		value=1;
		return;
	    else
		value=2;
		return;
	    end
	end

end


function dist = distanceTo(a,b)
    dist = (a-b) * (a-b);
end
    
	    
