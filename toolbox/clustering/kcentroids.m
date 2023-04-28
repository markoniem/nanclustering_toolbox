function [L, C, sumd] = kcentroids(X, k, varargin)
% Description:
% Performs K-means, K-medians or K-spatial-medians clustering. Three options
% are available for treating missing values. Weighted clustering is supported.  
%
% Inputs:
%               X - Input data set
%               k - Final number of clusters
%       'Weights' - Weights of observations. Default is a vector of ones.
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
%    'Replicates' - Selected number of repetitions. Default is 100.
%      'InitCrit' - Initialization criterion. Default is 'kmeans++'.
%                   Alternatives: 
%                   'random'   - random selection of initial prototypes 
%                   'kmeans++' - Uses K-means++ algorithm for selecting initial prototypes   
%        'Start' - Initial values of prototypes. Default is []
%'ShowProgression' - Indicator flag for presenting replicate number of clustering. Default value: false.
%   UseParallel' - Use parallelled computation. Default value: false.
%'CloseParallel' - Flag defines whether the parallel pool is closed. Default value: false.
%       'UseGPU' - Use GPU in computations. NOTE that the parrallelled GPU computation is not supported! 
%                  Default value: false.  
%
% Output:
%               L - Cluster labels for each observation
%               C - Cluster centroids
%            sumd - Sum of distances
%
pnames = {'weights' 'treatmissing' 'distance' 'replicates' 'initcrit' 'start' 'showprogression' 'useparallel' 'closeparallel' 'usegpu'};
dflts =  {ones(size(X,1),1) 'ads' 'euc' 100 'kmeans++' [] false false false false};
[weights, treatmissing, distance, replicates, initcrit, start, showprogression, useparallel, closeparallel, usegpu] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});
treatmissingNames = {'ads','pds','exp'};
distNames = {'sqe','euc','cit'};
initNames = {'random', 'kmeans++'};
weights = weights(:);
if ~isnumeric(weights) && size(X,1)~=length(weights)
    error('Invalid weights');
end
if isnumeric(treatmissing) || ~ismember(treatmissing,treatmissingNames)
    error('Invalid treatmissing');
end
if isnumeric(distance) || ~ismember(distance,distNames)
    error('Invalid distance');
end
if ~isnumeric(replicates)
    error('Invalid replicates');
end
if isnumeric(initcrit) || ~ismember(initcrit,initNames)
    error('Invalid initcrit');
end
if ~islogical(showprogression)
    error('Invalid showprogression');
end
if ~islogical(useparallel)
    error('Invalid useparallel');
end
if ~islogical(closeparallel)
    error('Invalid closeparallel');
end
if ~islogical(usegpu)
    error('Invalid usegpu');
end
if usegpu && useparallel
    error('Parallelised GPU is not supported');
end
if k > size(X,1)
    error('k is too high');
end
if sum(~any(isnan(X),2)) < k
    error('Not enough complete observations');
end
if (strcmp(treatmissing,'exp') && strcmp(distance,'cit'))
    error('Expected City block is not supported');
end
%
% Select distance computation functions
if ~any(isnan(X(:))) 
% No missing data. Use basic distance functions.    
        if strcmp(distance,'sqe')
            distFuncs = {@sqeucdistfuncp1; @sqeucdistfuncp2};
            centFun =  @meanfun;
        elseif strcmp(distance,'euc')
            distFuncs = {@eucdistfuncp1; @eucdistfuncp2};
            centFun = @spatialmedianfun;
        elseif strcmp(distance,'cit')
            distFuncs = {@citydistfuncp1; @citydistfuncp2};
            centFun = @medianfun;
        end    
else
% Available data strategy    
    if strcmp(treatmissing,'ads')
        if strcmp(distance,'sqe')
            distFuncs = {@sqeucdistfuncp1_ads; @sqeucdistfuncp2_ads};
            centFun =  @nanmeanfun;
        elseif strcmp(distance,'euc')
            distFuncs = {@eucdistfuncp1_ads; @eucdistfuncp2_ads};
            centFun =  @nanspatialmedianfun_ads;
        elseif strcmp(distance,'cit')
            distFuncs = {@citydistfuncp1_ads; @citydistfuncp2_ads};
            centFun =  @nanmedianfun;
        end  
% Partial distance strategy        
    elseif strcmp(treatmissing,'pds')
        if strcmp(distance,'sqe')
            distFuncs = {@sqeucdistfuncp1_pds; @sqeucdistfuncp2_pds};
            centFun =  @nanmeanfun;
        elseif strcmp(distance,'euc')
            distFuncs = {@eucdistfuncp1_pds; @eucdistfuncp2_pds};
            centFun =  @nanspatialmedianfun_pds;
        elseif strcmp(distance,'cit')
            distFuncs = {@citydistfuncp1_pds; @citydistfuncp2_pds};
            centFun =  @nanmedianfun;
        end        
% Expected distance strategy        
    elseif strcmp(treatmissing,'exp')
        if strcmp(distance,'sqe')
            distFuncs = {@sqeucdistfuncp1_exp; @sqeucdistfuncp2};
            centFun =  @nanmeanfun;
        elseif strcmp(distance,'euc')
            distFuncs = {@eucdistfuncp1_exp; @eucdistfuncp2_ads};
            centFun =  @nanspatialmedianfun_exp;
        end
    end        
end
%
if strcmp(treatmissing,'exp')
    % Compute parameters of conditional multivariate Gaussian distribution
    [Xi, sx, success] = ecmnmlefunc(X, 'Weights', weights);
    if ~success
        error('Error for computing parameters of conditional distribution');
    end
else
    Xi = zeros(size(X,1),1);
    sx = zeros(size(X,1),1);
end
%
if usegpu
    % Make sure that X, weights, and start are GPU arrays if usegpu is true  
    if ~isgpuarray(X)
        X = gpuArray(X);
    end
    if ~isgpuarray(weights)
        weights = gpuArray(weights);
    end    
    if ~isgpuarray(start)
        start = gpuArray(start);
    end
    if ~isgpuarray(Xi)
        Xi = gpuArray(Xi);
    end
    if ~isgpuarray(sx)
        sx = gpuArray(sx);
    end
end
%
% Use single CPU or GPU  
if ~useparallel
    L = [];
    C = [];
    sumd = Inf;
    sumdBest = Inf;
    for ii = 1:replicates
        if showprogression
            fprintf('Replicate: %d\n',ii);
        end
        [ClusterVal, success] = loopBody(X,k,weights,initcrit,start,distFuncs,centFun,usegpu,Xi,sx);
        if ~success
            if ~isempty(start)
                start(end,:) = [];
            end
        end
        L1 = ClusterVal{2};
        C1 = ClusterVal{3};
        sumd1 = ClusterVal{1};
        if sumd1 < sumdBest
            sumd = sumd1;
            L = L1;
            C = C1;
            sumdBest = sumd1;
        end
    end
else
    %
    % Use parallelled implementation
    L = cell(replicates,1);
    C = cell(replicates,1);
    sumd = cell(replicates,1);
    if isempty(gcp('nocreate'))
        parpool('local');
    end
    parfor ii = 1:replicates
        if showprogression
            fprintf('Replicate: %d\n',ii);
        end        
        [ClusterVal, ~] = loopBody(X,k,weights,initcrit,start,distFuncs,centFun,usegpu,Xi,sx);
        L{ii} = ClusterVal{2};
        C{ii} = ClusterVal{3};
        sumd{ii} = ClusterVal{1};
    end
    %
    % Find best value
    sumd = cell2mat(sumd);
    [~, idx] = min(sumd);
    sumd = sumd(idx);
    L = L{idx};
    C = C{idx};
end
%
% Close parallel pool 
if useparallel && closeparallel && ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

end

function [cellout, success] = loopBody(X,k,weights,initcrit,start,distFuncs,centFun,usegpu,Xi,sx)
% Description:
% Perform actual clustering process
%
% Inputs:
%             X - Input data set
%             k - Final number of clusters 
%       weights - Weights of observations
%      initcrit - Initialization criterion
%         start - Initial values of prototypes
%     distFuncs - Selected distance computation functions
%       centFun - Selected centroid computation function
%        usegpu - Use GPU in computations
%            Xi - Conditional mean imputed data set 
%                 (used only for computing expected distances)
%            sx - Variances of imputed data vectors 
%                 (used only for computing expected distances)  
%
% Output:
%    cellout{1} - Sum of distances
%    cellout{2} - Cluster labels for each observation
%    cellout{3} - Cluster centroids
%       success - Indicator flag for successfully finished algorithm  
%
success = 1;
distFuncp1 = distFuncs{1};
distFuncp2 = distFuncs{2};
cellout = cell(3,1);
%
% Initialization phase
weights1 = weights(~any(isnan(X),2));
X1 = X(~any(isnan(X),2),:);
if isempty(start), start = X1(randi(size(X1,1)),:); end
C = start;
if size(C,1) == 1
    if ~usegpu
        L = ones(size(X1,1),1);
    else
        L = ones(size(X1,1),1,"gpuArray");
    end
else
    L = distFuncp1(C,X1,Xi(~any(isnan(X),2),:),sx(~any(isnan(X),2),:));
end
if (strcmp(initcrit,'kmeans++'))
    for i = size(C,1)+1:k
        D = cumsum(weights1.*distFuncp2(C,X1,L));
        if D(end) == 0, success = 0; cellout{1} = Inf; return; end
        C(i,:) = X1(find(rand < D/D(end),1),:);
        L = distFuncp1(C,X1,Xi(~any(isnan(X),2),:),sx(~any(isnan(X),2),:));
    end
elseif (strcmp(initcrit,'random'))
    for i = size(C,1)+1:k
        C(i,:) = X1(randi(size(X1,1)),:);
    end
end
%
% Iterative refinement phase
L = distFuncp1(C,X,Xi,sx);
L1 = 0;
iter = 0;
while any(L ~= L1)
    L1 = L;
    for i = 1:k
        l = L==i;
        if sum(l) < 1
            success = 0;
        else
            C(i,:) = centFun(X(l,:),C(i,:),weights(l),Xi(l,:),sx(l,:));
        end
    end
    L = distFuncp1(C,X,Xi,sx);
    iter = iter + 1;
    if (iter > 100)
        warning('Maximum number of iterations occured in refinement phase!');
        break;
    end
end
if success
    D = weights.*distFuncp2(C,X,L);
    sumd = sum(D,'omitnan');
else
    sumd = Inf;
end
%
% Save results
cellout{1} = sumd;
cellout{2} = L;
cellout{3} = C;

end

function L = sqeucdistfuncp1(C, X, varargin) 
% Description: 
% Computes nearest centroids for each observation using squared Euclidean 
% distance.
%
% Inputs:
%        C - Matrix of cluster centroids
%        X - Input data set
%
% Output:
%        L - Cluster labels for each observation  
%
[~, L] = max(bsxfun(@minus,2*real(C*X'),dot(C,C,2)),[],1);
L = L(:);

end

function D = sqeucdistfuncp2(C, X, L, varargin) 
% Description: 
% Computes squared Euclidean distances between observations and centroids.
% 
% Inputs:
%        C - Matrix of cluster centroids
%        X - Input data set
%        L - Cluster labels for each observation
%
% Output:
%        D - Distances to nearest centroids
%
D = X-C(L,:);
D = sum(D.^2,2);

end

function L = sqeucdistfuncp1_ads(C, X, varargin) 
% Description: 
% Computes nearest centroids for each observation using squared Euclidean 
% distance. Missing values are treated using available data strategy.
%
% Inputs:
%        C - Matrix of cluster centroids
%        X - Input data set
%
% Output:
%        L - Cluster labels for each observation  
%
I = isnan(X);
X(I) = 0;
[~, L] = min((C.^2)*~I'-2*C*X',[],1);
L = L(:);

end

function D = sqeucdistfuncp2_ads(C, X, L) 
% Description: 
% Computes squared Euclidean distances between observations and centroids. 
% Missing values are treated using available data strategy.
% 
% Inputs:
%        C - Matrix of cluster centroids
%        X - Input data set
%        L - Cluster labels for each observation
%
% Output:
%        D - Distances to nearest centroids
%
D = X-C(L,:);
D = sum(D.^2,2,'omitnan');

end

function L = sqeucdistfuncp1_pds(C, X, varargin) 
% Description: 
% Computes nearest centroids for each observation using squared Euclidean 
% distance. Missing values are treated using partial distance strategy.
%
% Inputs:
%        C - Matrix of cluster centroids
%        X - Input data set
%
% Output:
%        L - Cluster labels for each observation  
%
I = isnan(X);
X(I) = 0;
[~, L] = min((C.^2)*~I'-2*C*X',[],1);
L = L(:);

end

function D = sqeucdistfuncp2_pds(C, X, L) 
% Description: 
% Computes squared Euclidean distances between observations and centroids. 
% Missing values are treated using partial distance strategy.
%
% Inputs:
%        C - Matrix of cluster centroids
%        X - Input data set
%        L - Cluster labels for each observation
%
% Output:
%        D - Distances to nearest centroids
%
D = X-C(L,:);
D = sum(D.^2,2,'omitnan');
D = D.*(size(X,2)./sum(~isnan(X),2));

end

function L = sqeucdistfuncp1_exp(C, ~, Xi, sx) 
% Description: 
% Computes nearest centroids for each observation using squared Euclidean
% distance. Missing values are treated using expected distance strategy.
%
% Inputs:
%         C - Matrix of cluster centroids
%        Xi - Inputed data set
%        sx - Variances of original data set  
%
% Output:
%         L - Cluster labels for each observation
%
I = ones(size(C));
I2 = ones(size(Xi));
D = ((Xi.^2)*I'-2*Xi*C'+I2*(C.^2)');
D = D + sum(sx,2);
[~, L] = min(D,[],2);
L = L(:);

end

function L = citydistfuncp1(C, X, varargin) 
% Description: 
% Computes nearest centroids for each observation using City block distance.
%
% Inputs:
%        C - Matrix of cluster centroids
%        X - Input data set
%
% Outputs:
%        L - Cluster labels for each observation
%
[nx, p] = size(X);
[nc, ~] = size(C);
D = zeros(nc,nx,'double');
for i = 1:nc
    dsq = zeros(nx,1,'double');
    for q = 1:p
        dsq1 = abs(X(:,q)-C(i,q));
        dsq = dsq + dsq1; 
    end
    D(i,:) = dsq;
end
[~, L] = min(D,[],1);
L = L(:);

end

function D = citydistfuncp2(C, X, L) 
% Description: 
% Computes City block distances between observations and centroids.
%
% Inputs:
%         C - Matrix of cluster centroids
%         X - Input data set
%         L - Cluster labels for each observation
%
% Output:
%         D - Distances to nearest centroids
%
D = X-C(L,:);
D = sum(abs(D),2);

end

function L = citydistfuncp1_ads(C, X, varargin) 
% Description: 
% Computes nearest centroids for each observation using City 
% block distance. Missing values are treated using available data strategy.
%
% Inputs:
%        C - Matrix of cluster centroids
%        X - Input data set
%
% Outputs:
%        D - Distances to nearest centroids
%        L - Cluster labels for each observation
%
[nx, p] = size(X);
[nc, ~] = size(C);
D = zeros(nc,nx,'double');
for i = 1:nc
    dsq = zeros(nx,1,'double');
    for q = 1:p
        dsq1 = abs(X(:,q)-C(i,q));
        dsq1(isnan(dsq1)) = 0;
        dsq = dsq + dsq1;
    end
    D(i,:) = dsq;
end
[~, L] = min(D,[],1);
L = L(:);

end


function D = citydistfuncp2_ads(C, X, L) 
% Description: 
% Computes City block distances between observations and centroids.
% Missing values are treated using available data strategy.
%
% Inputs:
%         C - Matrix of cluster centroids
%         X - Input data set
%         L - Cluster labels for each observation
%
% Output:
%         D - Distances to nearest centroids
%
D = X-C(L,:);
D = sum(abs(D),2,'omitnan');

end

function L = citydistfuncp1_pds(C, X, varargin) 
% Description: 
% Computes nearest centroids for each observation using City 
% block distance. Missing values are treated using partial distance strategy.
%
% Inputs:
%        C - Matrix of cluster centroids
%        X - Input data set
%
% Outputs:
%        D - Distances to nearest centroids
%        L - Cluster labels for each observation
%
[nx, p] = size(X);
[nc, ~] = size(C);
D = zeros(nc,nx,'double');
for i = 1:nc
    dsq = zeros(nx,1,'double');
    for q = 1:p
        dsq1 = abs(X(:,q)-C(i,q));
        dsq1(isnan(dsq1)) = 0;
        dsq = dsq + dsq1;
    end
    D(i,:) = dsq;
end
[~, L] = min(D,[],1);
L = L(:);

end

function D = citydistfuncp2_pds(C, X, L) 
% Description: 
% Compute City block distances between observations and centroids.
% Missing values are treated using partial distance strategy.
%
% Inputs:
%         C - Matrix of cluster centroids
%         X - Input data set
%         L - Cluster labels for each observation
%
% Output:
%         D - Distances to nearest centroids
%
D = X-C(L,:);
D = sum(abs(D),2,'omitnan');
D = D.*(size(X,2)./sum(~isnan(X),2));

end

function L = eucdistfuncp1(C, X, varargin) 
% Description: 
% Computes nearest centroids for each observation using Euclidean distance.
%
% Inputs:
%         C - Matrix of cluster centroids
%         X - Input data set
%
% Output:
%         L - Cluster labels for each observation  
%
[~, L] = max(bsxfun(@minus,2*real(C*X'),dot(C,C,2)),[],1);
L = L(:);

end

function D = eucdistfuncp2(C, X, L)
% Description: 
% Compute Euclidean distances between observations and centroids.
%
% Inputs:
%         C - Matrix of cluster centroids
%         X - Input data set
%         L - Cluster labels for each observation 
%
% Output:
%         D - Distances to nearest centroids 
%
D = X-C(L,:);
D = sqrt(sum(D.^2,2));

end

function L = eucdistfuncp1_ads(C, X, varargin) 
% Description: 
% Computes nearest centroids for each observation using Euclidean distance.
% Missing values are treated using available data strategy.
%
% Inputs:
%         C - Matrix of cluster centroids
%         X - Input data set
%
% Output:
%         L - Cluster labels for each observation  
%
I = isnan(X);
X(I) = 0;
[~, L] = min((C.^2)*~I'-2*C*X',[],1);
L = L(:);

end

function D = eucdistfuncp2_ads(C, X, L)
% Description: 
% Computes Euclidean distances between observations and centroids.
% Missing values are treated using available data strategy.
%
% Inputs:
%         C - Matrix of cluster centroids
%         X - Input data set
%         L - Cluster labels for each observation 
%
% Output:
%         D - Distances to nearest centroids 
%
D = X-C(L,:);
D = sqrt(sum(D.^2,2,'omitnan'));

end

function L = eucdistfuncp1_pds(C, X, varargin) 
% Description: 
% Computes nearest centroids for each observation using Euclidean distance.
% Missing values are treated using partial distance strategy.
%
% Inputs:
%         C - Matrix of cluster centroids
%         X - Input data set
%
% Output:
%         L - Cluster labels for each observation  
%
I = isnan(X);
X(I) = 0;
[~, L] = min((C.^2)*~I'-2*C*X',[],1);
L = L(:);

end

function D = eucdistfuncp2_pds(C, X, L)
% Description: 
% Computes Euclidean distances between observations and centroids.
% Missing values are treated using partial distance strategy.
%
% Inputs:
%         C - Matrix of cluster centroids
%         X - Input data set
%         L - Cluster labels for each observation 
%
% Output:
%         D - Distances to nearest centroids 
D = X-C(L,:);
D = sum(D.^2,2,'omitnan');
D = D.*(size(X,2)./sum(~isnan(X),2));
D = sqrt(D);

end

function L = eucdistfuncp1_exp(C, ~, Xi, sx)
% Description: 
% Computes nearest centroids for each observation using Euclidean 
% distance. Missing values are treated using expected distance strategy.
%
% Inputs:
%         C - Matrix of cluster centroids
%        Xi - Inputed data set
%        sx - Variances of missing values in original data X  
%
% Output:
%         L - Cluster labels for each observation
%
I = ones(size(C));
I2 = ones(size(Xi));
D = ((Xi.^2)*I'-2*Xi*C'+I2*(C.^2)');
D = D + sum(sx,2);
[~, L] = min(D,[],2);
L = L(:);

end


function C = meanfun(X, C, weights, varargin)
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
C(1,:) = sum(weights.*X,1) ./ sum(weights);

end

function C = nanmeanfun(X, C, weights, varargin)
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

function C = medianfun(X, ~, weights, varargin)
% Description: 
% Computes median value of input data set. 
%
% Input:
%         X - Input data set 
%   weights - Weights of input data    
%
% Output:
%         C - Median value of data set
%
n = size(X,2);
wMat = weights.*ones(size(X));
nn = sum(wMat,1);
for j = 1:n
    [X(:,j), idx] = sort(X(:,j));
    wMat(:,j) = wMat(idx,j); 
end
wMat = wMat./sum(wMat,1);
wMat(isinf(wMat)) = 0;
C = nan(1,n);
for j = 1:n
    jj = find(cumsum(wMat(:,j),1)>=0.5,1);
    cumsumW = cumsum(wMat(1:jj,j));
    cumsumW = round(gather(cumsumW(end)),10);
    if isempty(jj)
        C(1,j) = NaN;
    elseif mod(nn(j),2)==0 && cumsumW==0.5
        C(1,j) = 0.5*(X(jj,j)+X(jj+1,j));
    else
        C(1,j) = X(jj,j);
    end
end

end

function C = nanmedianfun(X, ~, weights, varargin)
% Description: 
% Computes median value of input data set. 
%
% Input:
%         X - Input data set 
%   weights - Weights of input data    
%
% Output:
%         C - Median value of data set
%
n = size(X,2);
wMat = weights.*~isnan(X);
nn = sum(wMat,1);
for j = 1:n
    [X(:,j), idx] = sort(X(:,j));
    wMat(:,j) = wMat(idx,j); 
end
wMat = wMat./sum(wMat,1);
wMat(isinf(wMat)) = 0;
C = nan(1,n);
for j = 1:n
    jj = find(cumsum(wMat(:,j),1)>=0.5,1);
    cumsumW = cumsum(wMat(1:jj,j));
    cumsumW = round(gather(cumsumW(end)),10);
    if isempty(jj)
        C(1,j) = NaN;
    elseif mod(nn(j),2)==0 && cumsumW==0.5
        C(1,j) = 0.5*(X(jj,j)+X(jj+1,j));
    else
        C(1,j) = X(jj,j);
    end
end

end

function C = spatialmedianfun(X, u, weights, varargin)
% Description:
% Computes spatial-median value of a cluster.
%
% Inputs:
%         X - Input data partition
%         u - Previous value of cluster centroid
%   weights - Weights of input data    
%
% Output:
%         C - Spatial-median value
%
max_iter = 100;
tol = 1e-5;
iters = 0;
w = 1.5;
while iters < max_iter
    iters = iters + 1;
    D = X-u;
    D = dot(D,D,2);
    a = weights./sqrt(D+sqrt(eps));
    ax = sum(a.*X,1);
    v = (1/sum(a)).*ax;
    u1 = u + w*(v-u);
    if norm(u1-u,inf) < tol
        break;
    end
    u = u1;
end
C = u1;

end

function C = nanspatialmedianfun_ads(X, u, weights, varargin)
% Description: 
% Computes spatial median value of input data set. 
% Uses available data strategy for treating missing values.
%
% Inputs: 
%         X - Input data set
%         u - Previous value of cluster centroid
%   weights - Weights of input data    
%
% Output: 
%         C - Spatial-median value of data set
%
max_iter = 100;
tol = 1e-5;
P = ~isnan(X);
iters = 0;
w = 1.5;
X(isnan(X)) = 0;
while iters < max_iter
    iters = iters + 1;
    D = P.*(X-u);
    a = weights./sqrt(sum(D.^2,2)+sqrt(eps));
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

function C = nanspatialmedianfun_pds(X, u, weights, varargin)
% Description: 
% Computes spatial median value of input data set. 
% Uses partial distance strategy for treating missing values.
%
% Inputs: 
%         X - Input data set
%         u - Previous value of cluster centroid
%   weights - Weights of input data    
%
% Output: 
%         C - Spatial-median value of data set
%
max_iter = 100;
tol = 1e-5;
P = ~isnan(X);
iters = 0;
w = 1.5;
Xorg = X;
X(isnan(X)) = 0;
while iters < max_iter
    iters = iters + 1;
    D = P.*(X-u);
    alpha = size(X,2)./sum(~isnan(Xorg),2);
    a = weights./sqrt(alpha.*sum(D.^2,2)+sqrt(eps));
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

function C = nanspatialmedianfun_exp(X, u, weights, Xi, sx)  
% Description: 
% Computes spatial median value of input data set. 
% Uses expected distance strategy for treating missing values.
%
% Inputs: 
%         X - Input data set
%         u - Previous value of cluster centroid
%   weights - Weights of input data    
%        Xi - Conditional mean imputed data set 
%        sx - Variances of imputed data vectors       
%
% Output: 
%         C - Spatial-median value of data set
%
max_iter = 100;
tol = 1e-5;
P = ~isnan(X);
iters = 0;
w = 1.5;
X(isnan(X)) = 0;
while iters < max_iter
    iters = iters + 1;
    D = Xi-u;
    a = weights./sqrt(sum(D.^2,2)+sum(sx,2)+sqrt(eps));
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


