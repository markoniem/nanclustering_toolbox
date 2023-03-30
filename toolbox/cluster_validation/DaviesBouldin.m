function DB = DaviesBouldin(X, C, L, weights, Xi, sx, varargin)
% Description: 
% Computes DaviesBouldin index value. Three options are available 
% for treating missing values. See documentation. 
%
% Inputs:
%               X - Input data set with missing values
%               C - Cluster centroids 
%               L - Cluster labels 
%         weights - Weights of observations. Default is a vector of ones.
%              Xi - Imputed data using conditional mean imputation
%              sx - Variance of data set
%  'TreatMissing' - Missing values treating strategy. Default is 'ads'.      
%                   Alternatives:
%                   'ads' - available data strategy
%                   'pds' - partial distance strategy
%                   'exp' - expected distance strategy                      
%      'Distance' - Selected distance metric. Default is 'euc'.
%                   Alternatives: 
%                   'sqe' - squared Euclidean distance
%                   'euc' - Euclidean distance
%                   'cit' - City block distance 
%                   Note: Expected City block distance is not supported.
%
% Output:
%        CH - Value of Calinski-Harabasz index 
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
        nanmatrixdist = @nanmatrixdist_ads;
        nanpdistfunc = @nanpdistfunc_ads;
    elseif strcmp(treatmissing,'pds')
        nanmatrixdist = @nanmatrixdist_pds; 
        nanpdistfunc = @nanpdistfunc_pds;
    end    
    clusts = unique(L);
    num = length(clusts);
    if num == 1, DB = Inf; return; end
    aveWithinD = zeros(num,1);
    for i = 1:num
        members = (L == clusts(i));
        aveWithinD(i) = nanmean(nanmatrixdist(X(members,:),C(i,:), ...
            distance),weights(members));
    end
    interD = nanpdistfunc(C,distance);
    R = zeros(num);
    cnt = 1;
    for i = 1:num
        for j = i+1:num
            R(i,j) = (aveWithinD(i)+aveWithinD(j)) / interD(cnt);
            cnt = cnt + 1;
        end
    end
    R = R + R';
    RI = max(R,[],1);
    DB = nanmean(RI(:),ones(length(RI(:)),1));
else
    clusts = unique(L);
    num = length(clusts);
    if num == 1, DB = Inf; return; end
    aveWithinD = zeros(num,1);
    for i = 1:num
        members = (L == clusts(i));
        aveWithinD(i) = nanmean(nanmatrixdist_exp(Xi(members,:), ...
            repmat(C(i,:),sum(members),1),distance,sx(members,:)),weights(members));
    end
    interD = nanpdistfunc_exp(C,distance);
    R = zeros(num);
    cnt = 1;
    for i = 1:num
        for j = i+1:num
            R(i,j) = (aveWithinD(i)+aveWithinD(j)) / interD(cnt);
            cnt = cnt + 1;
        end
    end
    R = R + R';
    RI = max(R,[],1);
    DB = nanmean(RI(:),ones(length(RI(:)),1));      
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

function D = nanmatrixdist_ads(X, Y, distance)
% Description: 
% Computes distances between two equal sized matrices using 
% available data strategy for treating missing values. 
%
% Inputs:
%         X - First data set
%         Y - Second data set
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance
%             'euc' - Euclidean distance 
%             'cit' - City block distance 
%             'sqcit' - squared City block distance 
%
% Output:
%         D - Distance matrix
switch distance
    case 'sqe'
        D = bsxfun(@minus,X,Y);       
        D = nansum(D.^2,2);
    case 'euc'
        D = bsxfun(@minus,X,Y);       
        D = sqrt(nansum(D.^2,2));
    case 'cit'
        D = bsxfun(@minus,X,Y);      
        D = nansum(abs(D),2);
    case 'sqcit'
        D = bsxfun(@minus,X,Y);
        D = nansum(abs(D).^2,2);
end
             
end

function D = nanmatrixdist_pds(X, Y, distance)
% Description: 
% Computes distances between two equal sized matrices using 
% partial distance strategy for treating missing values.
%
% Input:
%         X - First data set 
%         Y - Second data set
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance 
%             'euc' - Euclidean distance 
%             'cit' - City block distance  
%            'scit' - squared City block distance  
%
% Output:
%         D - Distance matrix
switch distance    
    case 'sqe'
        D = bsxfun(@minus,X,Y);
        m = sum(~isnan(D),2);        
        D = nansum(D.^2,2);
        D = D.*(size(X,2)./m);
    case 'euc'
        D = bsxfun(@minus,X,Y);
        m = sum(~isnan(D),2);        
        D = nansum(D.^2,2);
        D = sqrt(D.*(size(X,2)./m));
    case 'cit'
        D = bsxfun(@minus,X,Y);
        m = sum(~isnan(D),2);        
        D = nansum(abs(D),2);
        D = D.*(size(X,2)./m);
    case 'sqcit'
        D = bsxfun(@minus,X,Y);
        m = sum(~isnan(D),2);        
        D = nansum(abs(D),2);
        D = (D.*(size(X,2)./m)).^2;        
end
             
end

function D = nanmatrixdist_exp(X, Y, distance, sx, sy)
% Description: 
% Computes distances between two equal sized matrices using expected 
% distances for treating missing values. 
%
% Input:
%          X - First imputed data set
%          Y - Second imputed data set  
%   distance - Selected distance metric 
%              Alternatives: 
%              'euc' - Euclidean distance 
%              'sqe' - squared Euclidean distance
%         sx - Variance of X
%         sy - Variance of Y
%
% Output:
%         D - Distance matrix 
%
if nargin == 3
    sx = zeros(size(X));
    sy = zeros(size(Y));
end
if nargin == 4
    sy = zeros(size(Y));
end
switch distance
    case 'euc'    
        % Compute omega
        s = sum(sx,2) + sum(sy,2);
        omega = bsxfun(@minus,X,Y);       
        omega = sum(omega.^2,2);
        omega = omega + s;
        % Compute variance
        Ex = X;
        Ex2 = X.^2 + sx;
        Ex3 = X.^3 + 3*X.*sx;
        Ex4 = X.^4 + 6*(X.^2).*sx + 3*sx.^2;
        Ey = Y;
        Ey2 = Y.^2 + sy;
        Ey3 = Y.^3 + 3*Y.*sy;
        Ey4 = Y.^4 + 6*(Y.^2).*sy + 3*sy.^2;        
        var = sum(bsxfun(@plus,Ex4,Ey4) - 4*bsxfun(@times,Ex3,Ey) - ...
                    4*bsxfun(@times,Ex,Ey3) + 6*bsxfun(@times,Ex2,Ey2),2) - ...
                        sum((bsxfun(@plus,Ex2-2*bsxfun(@times,Ex,Ey),Ey2)).^2,2);
        var(var<0.0000001) = 0;
        % Compute EED
        m = (omega.^2)./var;
        D = exp(gammaln(m + 0.5) - gammaln(m));
        D = D.*((omega./m).^(0.5));
        ind = isnan(D);
        D(ind) = sqrt(omega(ind));
    case 'sqe'
        % Compute ESD
        sx = sum(sx,2); 
        sy = sum(sy,2);
        s = sx + sy;  
        D = bsxfun(@minus,X,Y);       
        D = sum(D.^2,2);
        D = D + s;            
end
             
end

function C = nanmean(X, weights, varargin)
% Description: 
% Computes mean value of input data set. 
%
% Input:
%         X - Input data set
%   weights - Weights of input data set  
%
% Output:
%         C - Mean value of data set  
%
nan = isnan(X);
X(nan) = 0;
C(1,:) = sum(weights.*X,1) ./ sum(weights.*~nan,1);

end

