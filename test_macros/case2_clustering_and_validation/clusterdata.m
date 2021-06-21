% Description: 
% Perform clustering and saves clustering results to file. The results includes 
% cluster centroids and cluster labels for each observation. 
%
% Inputs (see selected values from 'parameters.m' file):
%        datasets - Input data set(s)
%    datamatrices - Data sets which consist of predefined numbers of missing
%                   values (loaded from file by default)
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
clear, clc, close all;
% addpath('../../benchmark_data/O-sets/');
% addpath('../../benchmark_data/Sim-sets/');
% addpath('../../benchmark_data/S-sets/');
% addpath('../../toolbox/preprocess');
load('../../benchmark_data/CSAI_data/datamatrices');
addpath('../../toolbox/kcentroids');
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
%
%datamatrices = cell(length(datasets),length(probmiss));
centers = cell(maxClusters-1,length(datasets),length(probmiss));
labels = cell(maxClusters-1,length(datasets),length(probmiss));
for i = 1:length(datasets)
    fprintf('Data set: %s\n',datasets{i})
    for j = 1:length(probmiss)
        %load(char(datasets(i)));
        fprintf('Missing values: %.2f %% \n',probmiss(j));  
        %X = normalizedata(X,'min-max');
        %Xm = genmissdata(X,probmiss(j));
        %datamatrices{i,j} = Xm;
        Xm = datamatrices{i,j};
        [centers(:,i,j), labels(:,i,j)] = iterative_kcentroids(Xm,maxClusters, ...
          clustMethod,distance,replicates,useprevCent,initcrit,pipelined,showProgression); 
    end
    fprintf('=================\n');
end
%save('datamatrices','datamatrices');
save('centers','centers');
save('labels','labels');

