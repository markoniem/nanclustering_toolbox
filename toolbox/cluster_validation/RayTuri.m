function RT = RayTuri(X, C, L, weights, Xi, sx, varargin)
% Description: 
% Computes Ray-Turi index value. Three options are available 
% for treating missing values. See documentation. 
%
% Inputs:
%         X - Input data set with missing values
%         C - Matrix of cluster centroids 
%         L - Cluster labels for each observation
%   weights - Weights of observations. Default is a vector of ones.
%        Xi - Imputed data using conditional mean imputation
%        sx - Variance of data set
%  distance - Selected distance metric 
%             Alternatives: 
%             'euc' - Euclidean distance 
%             'sqe' - squared Euclidean distance
%             'cit' - City block distance  
%
% Output:
%       RT - Value of Ray-Turi index    
%
pnames = {'treatmissing' 'distance'};
dflts =  {'ads' 'euc'};
[treatmissing, distance] = internal.stats.parseArgs(pnames, dflts, varargin{:});
treatmissingNames = {'ads','pds','exp'};
distNames = {'sqe','euc','cit'};
if isnumeric(treatmissing) || ~ismember(treatmissing,treatmissingNames)
    error('Invalid treatmissing');
end
if isnumeric(distance) || ~ismember(distance,distNames)
    error('Invalid distance');
end
if (strcmp(treatmissing,'exp') && strcmp(distance,'cit'))
    error('Expected City block is not supported');
end
%
if ~strcmp(treatmissing,'exp')
    if strcmp(treatmissing,'ads')
        nandistfuncp2 = @nandistfuncp2_ads;
        nanpdistfunc = @nanpdistfunc_ads;
    elseif strcmp(treatmissing,'pds')    
        nandistfuncp2 = @nandistfuncp2_pds;
        nanpdistfunc = @nanpdistfunc_pds;    
    end
    K = size(C,1);
    if K == 1, RT = Inf; return; end
    Intra = nansum(weights.*nandistfuncp2(C,X,L,distance));
    Inter = min(nanpdistfunc(C,distance));
    RT = Intra / Inter;
else
    K = size(C,1);
    if K == 1, RT = Inf; return; end
    Intra = nansum(weights.*nandistfuncp2_exp(C,Xi,L,distance,sx));
    Inter = min(nanpdistfunc_exp(C,distance));
    RT = Intra / Inter;
end

end

function D = nandistfuncp2_ads(C, X, L, distance) 
% Description: 
% Computes distances between observations with missing values and nearest 
% centroids. Uses available data strategy for treating missing values. 
%
% Inputs:
%         C - Cluster centroids
%         X - Input data set
%         L - Cluster labels 
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance
%             'euc' - Euclidean distance 
%             'cit' - City block distance 
%             'sqcit' - squared City block distance 
%
% Output:
%         D -  Distances to nearest centroids 
%
switch distance
    case 'sqe'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(D.^2,2);
    case 'euc'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(D.^2,2);
        D = sqrt(D);        
    case 'cit'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(abs(D),2);
    case 'sqcit'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(abs(D).^2,2);        
end
        
end

function D = nandistfuncp2_pds(C, X, L, distance) 
% Description: 
% Computes distances between observations with missing values and nearest 
% centroids. Uses partial distance strategy for treating missing values. 
%
% Input:
%         C - Cluster centroids
%         X - Input data set
%         L - Cluster labels
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance 
%             'euc' - Euclidean distance 
%             'cit' - City block distance  
%            'scit' - squared City block distance  
%
% Output:
%         D -  Distances to nearest centroids
%
switch distance
    case 'sqe'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(D.^2,2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
    case 'euc'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(D.^2,2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
        D = sqrt(D);
    case 'cit'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(abs(D),2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
    case 'sqcit'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(abs(D),2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
        D = D.^2;
end
        
end

function D = nandistfuncp2_exp(C, Xi, L, distance, sx) 
% Description: 
% Computes distances between observations with missing values and nearest 
% centroids. Uses expected distance strategy for treating missing values. 
%
% Inputs:
%        C - Cluster centroids 
%       Xi - Imputed data set using conditional mean 
%        L - Cluster labels 
% distance - Selected distance metric 
%            Alternatives: 
%            'euc' - Euclidean distance 
%            'sqe' - squared Euclidean distance
%       sx - Variance of input data set
%
% Output:
%        D - Distances to nearest centroids
%
if nargin == 4
    sx = zeros(size(Xi));
end
switch distance
    case 'euc'
        % Compute omega
        omega = bsxfun(@minus,Xi,C(L,:));
        omega = sum(omega.^2,2) + sum(sx,2);  
        % Compute variance
        Ex = Xi;
        Ex2 = Xi.^2 + sx;
        Ex3 = Xi.^3 + 3*Xi.*sx;
        Ex4 = Xi.^4 + 6*(Xi.^2).*sx + 3*sx.^2;
        Y = C(L,:);
        Ey = Y;
        Ey2 = Y.^2;
        Ey3 = Y.^3;
        Ey4 = Y.^4;
        var = sum(Ex4 + Ey4 - 4*Ex3.*Ey - 4*Ex.*Ey3 + 6*Ex2.*Ey2,2) - ...
            sum((Ex2 - 2*Ex.*Ey + Ey2).^2,2);
        var(var<0.0000001) = 0; 
        % Compute EED
        m = (omega.^2)./var;
        D = exp(gammaln(m+0.5)-gammaln(m));
        D = D.*((omega./m).^(0.5));
        ind = isnan(D);
        D(ind) = sqrt(omega(ind));
    case 'sqe'
        % Compute ESD
        sx = sum(sx,2);
        D = bsxfun(@minus,Xi,C(L,:));
        D = sum(D.^2,2) + sx;
end

end

function D = nanpdistfunc_ads(X, distance)
% Description: 
% Computes pairwise distances between all pairs of observations in input data set. 
% Missing values are treated using available data strategy.
%
% Inputs:
%         X - Input data set 
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance
%             'cit' - City block distance 
%             'euc' - Euclidean distance 
%
% Output:
%           D - Distances between observations
%
switch distance
    case 'sqe'
        n = size(X,1);
        D = nandistfunc_ads(X,X,'sqe');
        B = tril(ones(n,n),-1);
        D = D(B==1);
    case 'cit'
        n = size(X,1);
        D = nandistfunc_ads(X,X,'cit');
        B = tril(ones(n,n),-1);
        D = D(B==1);
    case 'euc'
        n = size(X,1);
        D = nandistfunc_ads(X,X,'sqe');
        B = tril(ones(n,n),-1);
        D = sqrt(D(B==1));        
end

end

function D = nanpdistfunc_pds(X, distance)
% Description: 
% Computes pairwise distances between all pairs of observations in input data set. 
% Missing values are treated using partial distance strategy.
%
% Inputs:
%         X - Input data set 
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance
%             'cit' - City block distance 
%             'euc' - Euclidean distance 
%
% Output:
%           D - Distances between observations
%
switch distance
    case 'sqe'
        n = size(X,1);
        D = nandistfunc_pds(X,X,'sqe');
        B = tril(ones(n,n),-1);
        D = D(B==1);
    case 'cit'
        n = size(X,1);
        D = nandistfunc_pds(X,X,'cit');
        B = tril(ones(n,n),-1);
        D = D(B==1);
    case 'euc'
        n = size(X,1);
        D = nandistfunc_pds(X,X,'sqe');
        B = tril(ones(n,n),-1);
        D = sqrt(D(B==1));        
end

end

function D = nanpdistfunc_exp(X, distance, sx)
% Description: 
% Computes pairwise distances between all pairs of observations in data 
% matrix. Uses expected distance strategy for treating missing values. 
%
% Input:
%         X - Input data matrix 
%  distance - Selected distance metric 
%             Alternatives: 
%             'euc' - Euclidean distance 
%             'sqe' - squared Euclidean distance
%       sx - Variance of matrix 
%
% Output:
%        D - Distances between observations
%
if nargin == 2
    sx = zeros(size(X));
end
switch distance
    case 'euc'
        % Compute omega
        [N, ~] = size(X);
        s = sum(sx,2);
        ss = zeros(N*N-N*(N+1)/2,1);
        idx = 1;
        for i = 1:length(s)-1
            for j = i+1:length(s)
                ss(idx) = s(i) + s(j);
                idx = idx + 1;
            end
        end
        omega = (pdist(X)').^2 + ss;
        % Compute variance
        Ex = X;
        Ex2 = X.^2 + sx;
        Ex3 = X.^3 + 3*X.*sx;
        Ex4 = X.^4 + 6*(X.^2).*sx + 3*sx.^2;
        var1 = zeros(N*N-N*(N+1)/2,1);
        var2 = zeros(N*N-N*(N+1)/2,1);
        idx = 1;
        for i = 1:length(s)-1
            for j = i+1:length(s)
                var1(idx) = sum(Ex4(i,:) + Ex4(j,:) - 4*Ex3(i,:).*Ex(j,:) - ... 
                    4*Ex(i,:).*Ex3(j,:) + 6*Ex2(i,:).*Ex2(j,:),2);
                var2(idx) = sum((Ex2(i,:) - 2*Ex(i,:).*Ex(j,:) + Ex2(j,:)).^2,2);
                idx = idx + 1;
            end
        end
        var = var1 - var2;
        var(var<0.0000001) = 0;
        % Compute EED
        m = (omega.^2)./var;
        D = exp(gammaln(m+0.5) - gammaln(m));
        D = D.*((omega./m).^(0.5));
        ind = isnan(D);
        D(ind) = sqrt(omega(ind));   
    case 'sqe'
        % Compute ESD
        [N, ~] = size(X);
        s = sum(sx,2);
        ss = zeros(N*N-N*(N+1)/2,1);
        idx = 1;
        for i = 1:length(s)-1
            for j = i+1:length(s)
                ss(idx) = s(i) + s(j);
                idx = idx + 1;
            end
        end
        D = (pdist(X)').^2 + ss;    
end

end

function D = nandistfunc_ads(X, Y, distance)
% Description: 
% Computes pairwise distances between all pairs of observations in two 
% input data sets. Missing values are treated using available data strategy.
%
% Inputs:
%         X - First data set
%         Y - Second data set
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance
%             'cit' - City block distance 
%             'euc' - Euclidean distance 
%
% Output:
%         D - Distance matrix 
%
switch distance
    case 'sqe'
        I = isnan(Y);
        Y(I) = 0;
        I2 = isnan(X);
        X(I2) = 0;
        D = ((X.^2)*~I'-2*X*Y'+~I2*(Y.^2)');
        D(D<eps) = 0;
    case 'cit'
        [nx, p] = size(Y);
        [nc, ~] = size(X);
        D = zeros(nc,nx,'double');
        for i = 1:nc
            dsq = zeros(nx,1,'double');
            for q = 1:p
                dsq1 = abs(Y(:,q)-X(i,q));
                dsq1(isnan(dsq1)) = 0;
                dsq = dsq + dsq1;
            end
            D(i,:) = dsq;
        end
    case 'euc'
        I = isnan(Y);
        Y(I) = 0;
        I2 = isnan(X);
        X(I2) = 0;
        D = ((X.^2)*~I'-2*X*Y'+~I2*(Y.^2)');
        D(D<eps) = 0;
        D = sqrt(D);
end

end

function D = nandistfunc_pds(X, Y, distance)
% Description: 
% Computes pairwise distances between all pairs of observations in two 
% input data sets. Missing values are treated using partial distance strategy.
%
% Inputs: 
%        X - First data set
%        Y - Second data set
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance      
%             'cit' - City block distance  
%             'euc' - Euclidean distance 
%             
% Output:
%        D - Distance matrix 
switch distance
    case 'sqe'
        p = size(X,2);
        I = isnan(Y);
        Y(I) = 0;
        I2 = isnan(X);
        X(I2) = 0;
        D = ((X.^2)*~I'-2*X*Y'+~I2*(Y.^2)');
        D(D<eps) = 0;
        D = D.*(p./(double(~I2)*double(~I')));
    case 'cit'
        [nx, p] = size(Y);
        [nc, ~] = size(X);
        D = zeros(nc,nx,'double');
        for i = 1:nc
            dsq = zeros(nx,1,'double');
            m = p*ones(nx,1);
            for q = 1:p
                dsq1 = abs(Y(:,q)-X(i,q));
                isn = isnan(dsq1);
                m(isn) = m(isn) - 1;
                dsq1(isn) = 0;
                dsq = dsq + dsq1;
            end
            D(i,:) = dsq.*(p./m);
        end
    case 'euc'
        p = size(X,2);
        I = isnan(Y);
        Y(I) = 0;
        I2 = isnan(X);
        X(I2) = 0;
        D = ((X.^2)*~I'-2*X*Y'+~I2*(Y.^2)');
        D(D<eps) = 0;
        D = D.*(p./(double(~I2)*double(~I')));
        D = sqrt(D);
end

end

