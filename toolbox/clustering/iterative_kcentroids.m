function [centers, labels] = iterative_kcentroids(X, varargin)
% Description: 
% Performs kcentroids clustering iteratively from 2 to MaxClusters. 
% Centroids are based on the selected clustering algorithm.
%
% Inputs:
%                  X - Input data set
%      'MaxClusters' - Maximum number of clusters. Default value: 20 
%          'Weights' - Weights of observations. Default is a vector of ones.
%     'TreatMissing' - Missing values treating strategy. Default is 'ads'.      
%                      Alternatives:
%                      'ads' - available data strategy
%                      'pds' - partial distance strategy
%                      'exp' - expected distance strategy    
%         'Distance' - Selected distance metric. Default is 'euc'.
%                      Alternatives: 
%                      'sqe' - squared Euclidean distance
%                      'euc' - Euclidean distance
%                      'cit' - City block distance 
%                      Note: Expected City block distance is not supported.
%       'Replicates' - Selected number of repetitions. Default is 100.          
%         'InitCrit' - Initialization criterion. Default is 'kmeans++'.
%                      Alternatives: 
%                      'random'   - random selection of initial prototypes 
%                      'kmeans++' - Uses K-means++ algorithm for selecting initial prototypes
%      'UsePrevCent' - Use previous centers in initialization. Default value: true.   
%  'ShowProgression' - Indicator flag for presenting progression of clustering. Default value: true.
%        'Pipelined' - Pipeline clustering results based on expected distances. Default value: false.       
%      'UseParallel' - Use parallelled computation. Default value: false.
%    'CloseParallel' - Flag defines whether the parallel pool is closed. Default value: true.
%           'UseGPU' - Use GPU in computations. NOTE that the parrallelled GPU computation is not supported! 
%                      Default value: false.     
%
% Outputs:
%            centers - Obtained cluster centers
%             labels - Obtained cluster labels
%
pnames = {'maxclusters' 'weights' 'treatmissing' 'distance' 'replicates' 'initcrit' 'useprevcent' 'showprogression' 'pipelined' 'useparallel' 'closeparallel' 'usegpu'};
dflts =  {20 ones(size(X,1),1) 'ads' 'euc' 100 'kmeans++' true true false false true false};
[maxclusters, weights, treatmissing, distance, replicates, initcrit, useprevcent, showprogression, pipelined, useparallel, closeparallel, usegpu] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});
if (~isnumeric(maxclusters) || maxclusters < 2) 
    error('Invalid maxclusters');
end
if ~islogical(useprevcent)
    error('Invalid useprevcent');
end
if ~islogical(showprogression)
    error('Invalid showprogression');
end
if ~islogical(pipelined)
    error('Invalid pipelined');
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
centers = cell(maxclusters-1,1);
labels = cell(maxclusters-1,1);
C = [];
%
fprintf('Performing clustering...\n');
%
% Open parallel pool if needed
if useparallel && isempty(gcp('nocreate'))
    parpool('local');
end
%
for k = 2:maxclusters
    if showprogression
        if k~= maxclusters, fprintf('k: %d, ',k); else, fprintf('k: %d\n',k); end 
        if mod(k,10) == 0 && k~= maxclusters, fprintf('\n'); end
    end
    if ~useprevcent, C = []; end
    [L, C, ~] = kcentroids(X,k,'Weights',weights,'TreatMissing',treatmissing,'Distance',distance,'Replicates',replicates, ...
                        'InitCrit',initcrit,'Start',C,'UseParallel',useparallel,'UseGPU',usegpu);
    centers{k-1} = C;
    labels{k-1} = L;
end
%
% Pipelined clustering:
if (strcmp(treatmissing,'exp') && pipelined)
    for k = 2:maxclusters
        [labels{k-1}, centers{k-1}, ~] = ...
            kcentroids(X,k,'Weights',weights,'TreatMissing','ads','Distance',distance,'Replicates',1, ...
                        'InitCrit',initcrit,'Start',centers{k-1},'UseParallel',useparallel,'UseGPU',usegpu);
    end
end
fprintf('Done! \n');
%
% Close parallel pool 
if useparallel && closeparallel && ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

end

