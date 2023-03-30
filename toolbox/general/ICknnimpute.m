function Ximp = ICknnimpute(Xm,varargin)
% Description: 
% Performs Incomplete-Case k-nearest neighbors imputation for the data set with 
% missing values. 
%
% Inputs:
%             Xm - Original data set with missing values
%            'K' - Number of nearest neighbors. Default is: 5.   
% 'TreatMissing' - Missing values treating strategy. Default is 'ads'.      
%                  Alternatives:
%                  'ads' - available data strategy
%                  'exp' - expected distance strategy  
%     'Distance' - Selected distance metric. Default is 'euc'.
%                  Alternatives: 
%                  'sqe' - squared Euclidean distance
%                  'euc' - Euclidean distance
%                  'cit' - City block distance 
%                
% Outputs:
%        Ximp - Imputed data set 
%
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
[N, n] = size(Xm);
Im = find(any(isnan(Xm),2));
Ximp = Xm;
if (strcmp(treatmissing,'exp'))
    [Xi, sx, success] = ecmnmlefunc(Xm);
    if ~success
        error('Error for computing parameters of conditional distribution');
    end
end
for i = 1:N
    if ~ismember(i,Im)
        continue;
    end
    Xmi = Xm(i,:);
    miss = find(isnan(Xmi));
    notmiss = find(~isnan(Xmi));
    if strcmp(treatmissing,'ads')
        D = nandistfunc(Xmi,Xm,distance);
    else
        D = nandistfunc_exp(Xi(i,:), Xi, sx(i,:), sx);
    end 
    [~, idx] = sort(D);
    neighbors = zeros(length(miss),K);
    Kcounter = zeros(length(miss),1);
    for j = 2:length(idx)
        Xj = Xm(idx(j),:);
        miss2 = find(isnan(Xj));
        if sum(ismember(miss2,notmiss)) > 0
            continue;
        else
            for k = 1:length(miss)
                if ismember(miss(k),miss2) || Kcounter(k) == K
                    continue;
                else
                    Kcounter(k) = Kcounter(k) + 1;
                    neighbors(k, Kcounter(k)) = idx(j);
                end
            end
        end
        if isempty(find(neighbors==0,1))
            break;
        end
    end
    for j = 1:length(miss)
        if Kcounter(j) > 0
            neighbors1 = neighbors(j,1:Kcounter(j));
            Xmi(miss(j)) = nancentroid(Xm(neighbors1,miss(j)),distance);
        end
    end
    Ximp(i,:) = Xmi;
end
if ~isempty(isnan(Ximp))
    for i = 1:n
        X1 = Ximp(:,i);
        X1(isnan(X1)) = nancentroid(Xm(:,i),distance);
        Ximp(:,i) = X1;
    end
end

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

