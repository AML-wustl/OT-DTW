function out = clustereval(a, b, method)
%CLUSTEREVAL Easy clustering evaluation in MATLAB.
%   CLUSTEREVAL works with various evaluation metrics.
%
%   Example:
%
%       X = rand(100, 2);
%       Z = linkage(X, 'average', 'euclidean');
%       a = cluster(Z, 'maxclust', 4);
%       b = kmeans(X, 4);
%       clustereval(a, b, 'ari') % adjusted Rand index

%   Copyright (c) 2015 Taehoon Lee
%   ri: the Rand Index
%   Rand, "Objective Criteria for the Evaluation of Clustering Methods", JASA, 1971.
%   mi: the Mirkin index
%   hi: the Hubert index
%   ari: adjusted Rand index
%   Hubert and Arabie, "Comparing partitions", Journal of Classification, 1985.
%   fowlkes: the Fowlkes-Mallows index
%   Fowlkes and Mallows, "A Method for Comparing Two Hierarchical Clustering", JASA, 1983.
%   chi: Pearson's chi-square test
%   Chernoff and Lehmann, "The Use of Maximum Likelihood Estimates in \chi^2 Tests for Goodness of Fit", AMS, 1954.
%   cramer: Cramer's coefficient
%   tchouproff: Tchouproff's coefficient
%   moc: the Measure of Concordance
%   Pfitzner et al., "Characterization and evaluation of similaritymeasures for pairs of clusterings", KIS, 2009.
%   nmi: Normalized Mutual Index
%   Strehl and Ghosh, "Cluster ensembles - a knowledge reuse framework for combining multiple partitions", JMLR, 2002.

if nargin < 2 || min(size(a)) > 1 || min(size(b)) > 1 || numel(a) ~= numel(b)
   error('First two arguments must be two vectors with the same length.');
end

method = lower(method);
    n = numel(a);
    I = max(a);
    J = max(b);
    C = zeros(I, J);
    for i = 1:I
        tmp = a==i;
        for j = 1:J
            C(i,j) = sum(tmp & b==j);   % form contingency matrix
        end
    end

    nis = sum(sum(C,2).^2);             % sum of squares of sums of rows
    njs = sum(sum(C,1).^2);             % sum of squares of sums of columns
    
    nc2 = nchoosek(n, 2);               % all possible numbers
    n2 = sum(C(:).^2);                  % sum over rows & columnns of nij^2
    n3 = (nis + njs) / 2;               % no. agreements; 
                                        % disagreements - t2 + sumij*0.5;
                                        
    ni = sum(C,2);
    ni = ni.*(ni-1)/2;
    niss = sum(ni);
    nj = sum(C,1);
    nj = nj.*(nj-1)/2;
    njss = sum(nj);

    A = nc2 + n2 - n3;
    D = n3 - n2;
    
    CC = C.^2;
    CC = bsxfun(@rdivide, CC, sum(CC, 1));
    CC = bsxfun(@rdivide, CC, sum(CC, 2));
    chi = ( sum(sum(CC, 1), 2) - 1 ) * n;
    
    C = C / n;
    
    switch method
    case {'randindex', 'ri'}
        out = A / nc2;
    case {'mirkinindex', 'mi'}
        out = D / nc2;
    case {'hubertindex', 'hi'}
        out = (A - D) / nc2;
    case {'adjustedrandindex', 'adjustedri', 'adjustri', 'ari'}
        nc = (n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n) / (2*(n-1));
        if nc2 == nc
            out = 0;
        else
            out = (A-nc)/(nc2-nc);
        end
    case {'chisquared', 'chisquare', 'chi'}
        out = chi;
    case 'cramer'
        out = sqrt( chi / n / min(I-1, J-1) );
    case 'tchouproff'
        out = sqrt( chi / n / sqrt((I-1)*(J-1)) );
    case 'moc'
        if numel(C) == 1
            out = 1;
        else
            out = chi / n / ( sqrt(I*J) - 1 );
        end
    case 'nmi'
        Pxy = nonzeros(C);
        Px = sum(C, 1); % original: mean(C,1)... strange 
        Py = sum(C, 2);
        Hxy = -dot(Pxy, log2(Pxy+eps));
        Hx = -dot(Px, log2(Px+eps));
        Hy = -dot(Py, log2(Py+eps));
        MI = Hx + Hy - Hxy;
        out = sqrt((MI/Hx)*(MI/Hy));
    case 'jaccard'
        out = (n2 - n)/(nis + njs - n2 - n);	 % Jaccard - Jain and Dubes 1988
    case 'fm'
        out = 0.5*(njs - n)/sqrt(niss*njss);    % FM - Fowlkes and Mallows 1983
    end
    
end