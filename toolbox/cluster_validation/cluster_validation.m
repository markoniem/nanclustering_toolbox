function [indices_values, indices_names] = cluster_validation(X, centers, labels, varargin)
% Description: 
% Computes cluster validation indices values with missing values.
%
% Inputs:
%              X - Input data set with missing values
%        centers - Cluster centers obtained by iterative_kcentroids function
%         labels - Cluster labels obtained by iterative_kcentroids function
%      'Weights' - Weights of observations. Default is a vector of ones.
% 'TreatMissing' - Missing values treating strategy. Default is 'ads'.      
%                  Alternatives:
%                  'ads' - available data strategy
%                  'pds' - partial distance strategy
%                  'exp' - expected distance strategy                      
%     'Distance' - Selected distance metric. Default is 'euc'.
%                  Alternatives: 
%                  'sqe' - squared Euclidean distance
%                  'euc' - Euclidean distance
%                  'cit' - City block distance 
%                  Note: Expected City block distance is not supported.
%   'ClusterInd' - Selected cluster validation indices
%                  Alternatives (multiple indices can be selected): 
%                  @CalinskiHarabasz - Calinski-Harabasz (Selected by default)
%                  @DaviesBouldin - Davies-Bouldin (Selected by default)
%                  @DaviesBouldin2 - Davies-Bouldin* (Selected by default)
%                  @GenDunn - Generalized Dunn (Selected by default)
%                  @kCE - kCE-index (Selected by default)
%                  @PBM - Pakhira-Bandyopadhyay-Maulik (Selected by default)
%                  @RayTuri - Ray-Turi (Selected by default)
%                  @Silhouette - Silhouette (Selected by default)
%                  @Silhouette2 - Silhouette*
%                  @WB - WB-index (Selected by default)
%                  @WemmertGancarski - Wemmert-Gancarski (Selected by default) 
%                  %
%                  Note: Used distance metrics can be separately specified for 
%                  the selected indices as follows: 
%                  [{@index1, distance1}; {@index2, distance2}; i.e.]
%                  'Distance' option will be used if distances were not specified 
%                  separately for the indices (see definition above). 
%  'UseParallel' - Use parallelled computation. Default value: false.
%'CloseParallel' - Flag defines whether the parallel pool is closed. Default value: true.
%       'UseGPU' - Use GPU in computations. NOTE that the parrallelled GPU computation is not supported! 
%                  Default value: false.   
%
% Output:
%
% indices_values - Values of cluster validation indices.
% indices_names - Names of selected indices 
%
pnames = {'weights' 'treatmissing' 'distance' 'clusterind' 'useparallel' 'closeparallel' 'usegpu'};
dflts =  {ones(size(X,1),1)  'ads' 'euc' {@CalinskiHarabasz; @DaviesBouldin; @DaviesBouldin2; ... 
    @GenDunn; @kCE; @PBM; @RayTuri; @Silhouette; @WB; @WemmertGancarski}, false, true, false};
[weights, treatmissing, distance, clusterind, useparallel, closeparallel, usegpu] = ... 
                                        internal.stats.parseArgs(pnames, dflts, varargin{:});
treatmissingNames = {'ads','pds','exp'};
distNames = {'sqe','euc','cit'};
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
if (strcmp(treatmissing,'exp') && strcmp(distance,'cit'))
    error('Expected City block is not supported');
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
%
indices_names = clusterind;
indices_values = zeros(length(clusterind), size(centers,1)+1);
I = ~all(isnan(X),2);
X = X(I,:);
fprintf('Performing cluster validation...\n');
%
% Open parallel pool if needed
if useparallel && isempty(gcp('nocreate'))
    parpool('local');
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
    if ~isgpuarray(centers{1,1})
        for i = 1:length(centers)
            centers{i,1} = gpuArray(centers{i,1});
        end
    end
    if ~isgpuarray(labels{1,1})
        for i = 1:length(labels)
            labels{i,1} = gpuArray(labels{i,1});
        end
    end    
end
%
% Use single CPU or GPU  
if ~useparallel
    for k = 2:size(centers,1)+1
        C = centers{k-1};
        L = labels{k-1};
        L(~I) = 0;
        L = L(L~=0);
        L = L(:);
        if (strcmp(treatmissing,'ads') || strcmp(treatmissing,'pds'))
            for r = 1:size(clusterind,1)
                clusterIndex = clusterind{r,1};
                indices_values(r,1) = Inf(1,1);
                if size(clusterind,2) > 1
                    indexDist = clusterind{r,2};
                    indices_values(r,k) = clusterIndex(X,C,L,weights,[],0,'TreatMissing',treatmissing,'Distance',indexDist);
                else
                    indices_values(r,k) = clusterIndex(X,C,L,weights,[],0,'TreatMissing',treatmissing,'Distance',distance);
                end
            end
        elseif (strcmp(treatmissing,'exp'))
            if ~usegpu
                Xi = zeros(size(X));
                sx = zeros(size(X));                
            else
                Xi = gpuArray(zeros(size(X)));
                sx = gpuArray(zeros(size(X)));
            end
            for i = 1:k
                try
                    [Xi(L==i,:), sx(L==i,:), ~] = ecmnmlefunc(X(L==i,:),'Weights',weights(L==i));
                catch
                    Xi(L==i,:) = X(L==i,:);
                    sx(L==i,:) = zeros(size(sx(L==i,:)));
                    for j = 1:size(Xi,2)
                        X1 = Xi(L==i,j);
                        X1(isnan(X1)) = median(X1,'omitnan');
                        Xi(L==i,j) = X1;
                    end
                end
            end
            for r = 1:size(clusterind,1)
                clusterIndex = clusterind{r,1};
                indices_values(r,1) = Inf(1,1);
                if size(clusterind,2) > 1
                    indexDist = clusterind{r,2};
                    indices_values(r,k) = clusterIndex(X,C,L,weights,Xi,sx,'TreatMissing',treatmissing,'Distance',indexDist);
                else
                    indices_values(r,k) = clusterIndex(X,C,L,weights,Xi,sx,'TreatMissing',treatmissing,'Distance',distance);
                end
            end
        end
    end
%
% Use parallelled implementation
else
    numOfIndices = length(clusterind);
    indices_values = Inf(length(clusterind),size(centers,1)+1);
    indices_values(:,1) = Inf(size(indices_values,1),1);
    parfor k = 2:size(centers,1)+1
        C = centers{k-1};
        L = labels{k-1};
        L(~I) = 0;
        L = L(L~=0);
        L = L(:);
        if (strcmp(treatmissing,'ads') || strcmp(treatmissing,'pds'))
            for r = 1:numOfIndices
                clusterIndex = clusterind{r,1};
                if size(clusterind,2) > 1
                    indexDist = clusterind{r,2};
                    indices_values(r,k) = clusterIndex(X,C,L,weights,[],0,'TreatMissing',treatmissing,'Distance',indexDist);
                else
                    indices_values(r,k) = clusterIndex(X,C,L,weights,[],0,'TreatMissing',treatmissing,'Distance',distance);
                end
            end
        elseif (strcmp(treatmissing,'exp'))
            Xi = zeros(size(X));
            sx = zeros(size(X));
            for i = 1:k
                try
                    [Xi(L==i,:), sx(L==i,:), ~] = ecmnmlefunc(X(L==i,:),'Weights',weights(L==i));
                catch
                    Xi(L==i,:) = X(L==i,:);
                    sx(L==i,:) = zeros(size(sx(L==i,:)));
                    for j = 1:size(Xi,2)
                        X1 = Xi(L==i,j);
                        X1(isnan(X1)) = median(X1,'omitnan');
                        Xi(L==i,j) = X1;
                    end
                end
            end
            for r = 1:numOfIndices
                clusterIndex = clusterind{r,1};
                if size(clusterind,2) > 1
                    indexDist = clusterind{r,2};
                    indices_values(r,k) = clusterIndex(X,C,L,weights,Xi,sx,'TreatMissing',treatmissing,'Distance',indexDist);
                else
                    indices_values(r,k) = clusterIndex(X,C,L,weights,Xi,sx,'TreatMissing',treatmissing,'Distance',distance);
                end
            end
        end
    end 
end
fprintf('Done! \n');
%
% Close parallel pool 
if useparallel && closeparallel && ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

end

