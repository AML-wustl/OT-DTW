function DM = DMComputation(DataSet, DistanceIndex)
    % Distance Measures into consideration
    % ED = 1
    % NCCc = 2
    % cDTW with defined window = 3
    % LOROT_Wasserstein = 4
%     Distances = [cellstr('ED'), 'SBD', 'DTW']; 

    % For cDTW
    [~, len] = size(DataSet);
    DTW_Window = len;
    lambda_1 = 1;
    lambda_2 = 1;
    delta = 1;
    sigma = 2;

    DM = [];

    switch DistanceIndex            
        case 1
            DM = ParallelDM_ED(DataSet);
        case 2
            DM = ParallelDM_SBD(DataSet);
        case 3
            DM = ParallelDM_DTW(DataSet,DTW_Window);
        case 4
            DM = ParallelDM_LRW(DataSet,lambda_1,lambda_2,delta,sigma);
    end

%     dlmwrite( strcat( 'DATASETS/',char(Datasets(i)),'/', char(Datasets(i)), '_', char(Distances(DistanceIndex)) ,'.distmatrix'), DM, 'delimiter', '\t');
end