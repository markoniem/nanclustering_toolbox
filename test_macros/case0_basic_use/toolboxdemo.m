% Description: 
% Demostrates the basic use of nanclustering_toolbox. Includes data 
% pre-processing, clustering, and cluster validation.     
%
clear, clc, close all;
%
% LOAD BENCHMARK DATA 
load fisheriris;
X = meas;
clear meas; clear species;
%X = X(:,1:2);
%
% DATA PREPROCESSING
addpath('../../toolbox/preprocess');
addpath('../../toolbox/distance_strategies');
% Generate missing values
%help genmissdata;
Xm = genmissdata(X, 0.1);
% Normalize data
%help normalizedata
Xnorm = normalizedata(Xm, 'min-max', [-1, 1]);
% Perform k-nearest neighbors imputation
%help knnimpute
%Ximp = knnimpute(Xnorm, 5, 'sqe', 'ads');
Ximp = knnimpute(Xnorm, 5, 'sqe', 'exp');
% Transform data set to more spherical symmetric
%help datasetmap
Xmapped = datasetmap(Ximp);
scatter_data(Xmapped);
% CLUSTERING
addpath('../../toolbox/kcentroids');
% Cluster data
%help kcentroids
[L, C, sumd] = kcentroids(Xnorm, 5, 100, 'euc', 'kmeans++', []);
%
% CLUSTER VALIDATION
addpath('../../toolbox/cluster_indices');
% Use iterative clustering function to obtain centroids and labels
%help iterative_kcentroids;
[centers, labels] = iterative_kcentroids(Xm, 10);
% Define the used cluster validation indices
%help cluster_validation
dist = 'euc';
indices = {@CalinskiHarabasz; @DaviesBouldin; @kCE; @PBM; @RayTuri; ... 
                @Silhouette; @WB; @WemmertGancarski};
% Optional way to define the indices
% indices = [{@CalinskiHarabasz, 'sqe'}; {@DaviesBouldin, 'euc'}; {@kCE, 'sqe'}; ... 
%             {@PBM, 'euc'}; {@RayTuri, 'sqe'}; {@Silhouette, 'euc'}; {@WB, 'sqe'}; ... 
%             {@WemmertGancarski, 'euc'}];
% Perform actual cluster validation
indices_values = cluster_validation(Xnorm, centers, labels, dist, indices, 'ads');
%indices_values = cluster_validation(Xnorm, centers, labels, 'sqe', indices, 'exp');
% Visualize the results
plot_indices(indices, indices_values);
%
% Restore default path
restoredefaultpath;

