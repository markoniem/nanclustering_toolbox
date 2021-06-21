% Description: 
% Produce visualizations based on clustering and cluster validation 
% indices results  
%
% Inputs (see used values from 'parameters.m' file): 
%
%         datamatrices - Clustered data set(s)
%              centers - Obtained cluster centroids 
%               labels - Obtained cluster labels for each observation     
%             datasets - Names of data set(s)
%       indices_values - Calculated index values for each number of clusters. 
%             probmiss - Probabilities of missing values generated to data set(s)
%           clusterInd - Used cluster validation indices   
%
clear, clc, close all;
load('../../benchmark_data/CSAI_data/datamatrices');
load('centers');
load('labels');
load('indices_values');
addpath('../../toolbox/kcentroids');
addpath('../../toolbox/cluster_indices');
params = params();
datasets = params.datasets;  
probmiss = params.probmiss;
clusterInd = params.clusterInd;
%
% Visualize clustering result of O2000 data which consists of 10 % 
% missing values. Use K = 5 as selected number of clusters. 
dataset = 'O200';
missing_values = 0.00;
K = 5;
idx1 = find(matches(datasets,dataset));
idx2 = find(probmiss==missing_values);
X = datamatrices{idx1,idx2};
C = centers{K-1,idx1,idx2};
L = labels{K-1,idx1,idx2};
scatter_results(X,C,L);
%
% Plot cluster validatin curves for same data set.
indices_values = indices_values(:,:,idx1,idx2);
plot_indices(clusterInd,indices_values);

