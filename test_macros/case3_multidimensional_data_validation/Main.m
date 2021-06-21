% Description: 
% Perform clustering and cluster validation. Saves clustering and cluster 
% validation index results to files. The results includes cluster centroids, 
% cluster labels for each observation, and cluster validation indices measures. 
%
% Inputs (see current values from 'parameters.m' file):
%        datasets - Input data set(s)
%        probmiss - Probabilities of missing values generated to data set(s)
%     maxClusters - Maximum number of clusters
%     clustMethod - Selected clustering algorithm. 
%                   Alternatives: 
%                   @kcentroids - Available data based clustering 
%                   @kcentroids_partial - Partial distances based clustering 
%                   @kcentroids_expected - Expected distances based clustering
%        distance - Selected distance metric. 
%                   Alternatives: 
%                   'sqe' - squared Euclidean distance
%                   'euc' - Euclidean distance
%                   'cit' - City block distance 
%                   Note: @kcentroids_expected does not support 'cit' option 
%      replicates - Number of restarts of clustering
%     useprevCent - Boolean value that defines if previously saved
%                   centroids are used as initial values for clustering
%        initcrit - Initialization criterion
%                   Alternatives: 
%                   'random' - Random selection of initial points
%                   'kmeans++' - Kmeans++ based selection of initial points
%       pipelined - Pipeline expected distances based clustering results 
%                   Note: @kcentroids_expected only uses pipeline option   
% showProgression - Boolean value which indicate if progression of clustering   
%                   is presented. Default value: true
%     distFuncInd - Used distance estimation strategy for indices 
%                   Alternatives:
%                   'ads' - Available data 
%                   'pds' - Partial distance 
%                   'exp' - Expected distance (ESD/EED)
%      clusterInd - Selected cluster validation indices
%                   Alternatives (multiple indices can be selected): 
%                   @CalinskiHarabasz - Calinski-Harabasz 
%                   @DaviesBouldin - Davies-Bouldin 
%                   @DaviesBouldin2 - Davies-Bouldin*
%                   @GenDunn - Generalized Dunn
%                   @kCE - kCE-index
%                   @PBM - Pakhira-Bandyopadhyay-Maulik
%                   @RayTuri - Ray-Turi
%                   @Silhouette - Silhouette
%                   @Silhouette2 - Silhouette*
%                   @WB - WB-index
%                   @WemmertGancarski - Wemmert-Gancarski
%                   Note: Used distance metric(s) for the selected indices can 
%                   be optionally separately specified. 
%                   Use then the following format typing: 
%                   [{@index1, distance1}; {@index2, distance2}; ...], e.g.,
%                   [{@CalinskiHarabasz, 'sqe'}; {@DaviesBouldin, 'euc'}]
%                   In general case, if specific distances were not given, 
%                   the 'distance' parameter will be used (see definition above). 
%
% Outputs:
%         centers - Cluster centers. Centers will be saved to three 
%                   dimensional cell array. First dimension is for different 
%                   cluster numbers, second dimension is for different data 
%                   sets, and third dimension is for different percentages 
%                   of missing values in data.
%          labels - Cluster labels for each observation. Labels will be saved 
%                   to three dimensional cell array. First dimension is for 
%                   different cluster numbers, second dimension is for 
%                   different data sets, and third dimension is 
%                   for different percentages of missing values in data.
%  indices_values - Calculated index values for each number of clusters. 
%                   Indices will be saved to four dimensional matrix 
%                   such as first dimension is for cluster validation indices, 
%                   second dimension is for different number of clusters 
%                   (starting from 2), third dimension correspond different 
%                   data sets, and fourth dimension is for different number 
%                   of missing values.
% DyymmddThhm.xlsx - Recommended number of clusters by cluster validation 
%                    indices saved to .xlsx file format. File is named based 
%                    on execution time of macro.  
%
clear, clc, close all; 
addpath('../../benchmark_data/M_spheres')
addpath('../../toolbox/preprocess');
addpath('../../toolbox/kcentroids');
addpath('../../toolbox/cluster_indices');
params = params();
datasets = params.datasets;
maxClusters = params.maxClusters;
clustMethod = params.clustMethod;
distance = params.distance;
replicates = params.replicates;
useprevCent = params.useprevCent;
initcrit = params.initcrit;
pipelined = params.pipelined;
probmiss = params.probmiss;
showProgression = params.showProgression;
distFuncInd = params.distFuncInd;
clusterInd = params.clusterInd;
%
datamatrices = cell(length(datasets),length(probmiss));
centers = cell(maxClusters-1,length(datasets),length(probmiss));
labels = cell(maxClusters-1,length(datasets),length(probmiss));
indices_values = zeros(length(clusterInd),maxClusters, ... 
                    length(datasets),length(probmiss)); 
for i = 1:length(datasets)
    load(datasets{i});
    clear C;
    clear labels;
    fprintf('Data set: %s\n',datasets{i})
    for j = 1:length(probmiss)
        fprintf('Missing values: %.2f %% \n',probmiss(j));  
        X = normalizedata(X,'min-max');
        Xm = genmissdata(X,probmiss(j));    
        % At least 0.5 % of complete observations are required 
        p = uint64(0.005*size(Xm,1));
        I1 = ~any(isnan(Xm),2);
        if sum(I1) < p
            I1 = randperm(size(X,1));
            Xm(I1(1:30),:) = X(I1(1:p),:);
        end
        datamatrices{i,j} = Xm;
        % Data set mapping function
        %Xm = datasetmap(Xm);
        [centers(:,i,j), labels(:,i,j)] = iterative_kcentroids(Xm,maxClusters, ...
          clustMethod,distance,replicates,useprevCent,initcrit,pipelined,showProgression); 
        indices_values(:,:,i,j) = cluster_validation(Xm,centers(:,i,j), ... 
                 labels(:,i,j),distance,clusterInd,distFuncInd);
    end
    fprintf('=================\n');
end
save('datamatrices','datamatrices');
save('centers','centers');
save('labels','labels');
save('indices_values','indices_values');
%
% Save results to .xlsx file
%
[~, indices_results]= min(indices_values,[],2);
results = zeros(size(clusterInd,1),length(datasets),length(probmiss));
for i = 1:length(datasets), for j = 1:length(probmiss), results(:,i,j) = ...
            indices_results(:,:,i,j); end, end
for i = 1:length(datasets), datasets{i} = sprintf('Data set: %s',datasets{i}); end
clusterIndchar = cell(size(clusterInd,1),1);
for i = 1:size(clusterInd,1), clusterIndchar{i} = sprintf('%s',char(clusterInd{i,1})); end
for i = 1:size(results,3)    
    R = results(:,:,i);
    T = table(clusterIndchar,R);
    resultsTable = splitvars(T);
    resultsTable.Properties.VariableNames = [' ', datasets'];
    writetable(resultsTable,sprintf('D%sT%s.xlsx', ...
                datestr(now,'yymmdd'), datestr(now,'HHMM')), 'Sheet', ... 
                        sprintf('Missing values %2d%%',100*probmiss(i)));
end
restoredefaultpath;