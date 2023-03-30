function [D, success] = kNNI(Xm,varargin)
% Description: 
% Performs k-nearest neighbors imputation for data set with 
% missing values and returns pairwise distances between observatorions.  
%
% Inputs:
%              Xm - Original data set with missing values
%             'K' - Number of nearest neighbors. Imputation is performed 
%                   in complete manner if K = 1. Otherwise, centroid values of 
%                   k-nearest neighbors will be used. Default is: 5.  
%   'TreatMissing' - Missing values treating strategy. Default is 'ads'.      
%                   Alternatives:
%                   'ads' - available data strategy
%                   'exp' - expected distance strategy  
%       'Distance' - Selected distance metric. Default is 'euc'.
%                   Alternatives: 
%                   'sqe' - squared Euclidean distance
%                   'euc' - Euclidean distance
%                   'cit' - City block distance 
%                   Note: Expected City block distance is not supported.
%                
% Outputs:
%        Ximp - Imputed data set 
%     success - Indicator flag for successfully finished algorithm  
%
success = 1;
pnames = {'K' 'treatmissing' 'distance'};
dflts =  {5 'ads' 'euc'};
[K, treatmissing, distance] = internal.stats.parseArgs(pnames, dflts, varargin{:});
treatmissingNames = {'ads','exp'};
distNames = {'sqe','euc','cit'};
if isnumeric(treatmissing) || ~ismember(treatmissing,treatmissingNames)
    error('Invalid treatmissing');
end
if isnumeric(distance) || ~ismember(distance,distNames)
    error('Invalid distance');
end
if (~isnumeric(K) || K < 1)
    error('Invalid K value');
end
if (strcmp(treatmissing,'exp') && strcmp(distance,'cit'))
    error('Expected City block is not supported');
end
%
I = ~all(isnan(Xm),2);
Xm = Xm(I,:);
Ximp = Xm;
Im = find(any(isnan(Xm),2));
Ic = find(~any(isnan(Xm),2)); 
if (strcmp(treatmissing,'exp'))
    [Xi, sx, success] = ecmnmlefunc(Xm);
end
switch K
    % Nearest complete observation
    case 1
        for i = 1:length(Im)
            Xmi = Xm(Im(i),:);
            if strcmp(treatmissing,'ads')
                D = nandistfunc(Xmi,Xm(Ic,:),distance); 
            else
                D = nandistfunc_exp(Xi(Im(i),:), Xi(Ic,:), sx(Im(i),:), sx(Ic,:));
            end
            [~, idx] = min(D);
            nearest = Xm(Ic(idx),:);
            ind = find(isnan(Xmi));
            Xmi(ind) = nearest(ind);
            Ximp(Im(i),:) = Xmi;
        end
    % Centroid of K nearest observations
    otherwise
        for i = 1:length(Im)
            Xmi = Xm(Im(i),:);
            if strcmp(treatmissing,'ads')
                D = nandistfunc(Xmi,Xm,distance);
            else
                D = nandistfunc_exp(Xi(Im(i),:), Xi, sx(Im(i),:), sx);
            end
            [~, idx] = sort(D);
            cnt = 0;
            % In the case data vector consists missing values after 
            % imputation, K will be increased one by one until vector is 
            % fully complete. 
            while sum(isnan(Xmi)) > 0
                knearests = Xm(idx(2:K+1+cnt),:);
                centroid = nancentroid(knearests,distance);
                ind = find(isnan(Xmi));
                Xmi(ind) = centroid(ind);
                Ximp(Im(i),:) = Xmi;
                cnt = cnt + 1;
            end
        end
end
D = pdist(Ximp)';

end

function D = nandistfunc(X, Y, distance)
% Description: 
% Estimates distances between observations with missing values in two data 
% matrices. Available data strategy is used for treating missing values.
%
% Function call:
%        D = nandistfunc(X, Y, distance)
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

function D = nandistfunc_exp(X, Y, sx, sy)
% Description: 
% Computes distances between two matrices. Uses expected distance 
% estimation for treating missing values. 
%
% Inputs: 
%         X - Imputed data set using conditional mean
%         Y - Second imputed data set using conditional mean
%        sx - Variance of data set X
%        sy - Variance of data set Y
%           
% Output:
%       D - Distance matrix
%
if nargin == 3
    sx = zeros(size(X));
    sy = zeros(size(Y));
end
if nargin == 4
    sy = zeros(size(Y));
end
% Compute ESD
D = pdist2(X,Y).^2;
sx = sum(sx,2);
sy = sum(sy,2);
D = bsxfun(@plus,D,sx);
D = bsxfun(@plus,D,sy');
    
end

function C = nancentroid(X, distance)
% Description: 
% Computes centroid value of input data set. Algorithm uses available data 
% strategy for treating missing values.   
% 
% Inputs: 
%         X - Input data set 
%  distance - Selected distance metric  
%             Alternatives:
%             'euc' - Euclidean distance
%             'sqe' - squared Euclidean distance
%             'cit' - City block distance
%
% Output: 
%         C - Centroid value of data set 
%
switch distance
    case 'euc'
        % Spatial median value of data
        C = nanspatialmedianfunc(X);
    case 'sqe'
        nan = isnan(X);
        X(nan) = 0;
        % Mean value of data 
        C(1,:) = sum(X,1) ./ sum(~nan,1);  
    case 'cit'
        Xsorted = sort(X,1);
        count = sum(~isnan(X),1);
        nn = floor(0.5*count);
        n = size(X,2);
        C = NaN(1,n);
        % Median value of data 
        for j = 1:n
            if count(j) == 0
                C(1,j) = NaN;
            elseif mod(count(j),2) == 0
                C(1,j) = 0.5*(Xsorted(nn(j),j)+Xsorted(nn(j)+1,j));
            else
                C(1,j) = Xsorted(nn(j)+1,j);
            end
        end
end

end

function C = nanspatialmedianfunc(X)
% Description: 
% Computes spatial median value of input data set. 
% Uses available data strategy for treating missing values.
%
% Inputs: 
%         X - Input data set
%
% Output: 
%         C - Spatial-median value of data set
%
u = zeros(1,size(X,2));
for i = 1:size(X,2)
    I = ~isnan(X(:,i));
    u(i) = median(X(I,i));
end
max_iter = 100;
tol = 1e-5;
P = ~isnan(X);
iters = 0;
w = 1.5;
X(isnan(X)) = 0;
while iters < max_iter
    iters = iters + 1;
    D = P.*(X-u);
    a = 1./sqrt(sum(D.^2,2)+sqrt(eps));
    a = P.*a;
    ax = sum(a.*X,1);
    v = (1./sum(a)).*ax;
    u1 = u + w*(v-u);
    if norm(u1-u,inf) < tol
        break;
    end
    u = u1;
end
C = u1;

end

