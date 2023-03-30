% Description: 
% Demostrates the basic use of nanclustering_toolbox. Includes data 
% pre-processing, clustering, and cluster validation steps.     
%
clear, clc, close all;
addpath('../toolbox/general');
rng(1); %for reproducibility
%
% LOAD BENCHMARK DATA 
addpath('../datasets');
addpath('../toolbox/general');
load('O2000.mat');
%
% DATA PREPROCESSING
% Generate missing values
%help genmissdata;
[Xm, I] = genmissdata(X, 0.1);
Xm = Xm(I,:);
% Normalize data
%help normalizedata
Xnorm = normalizedata(Xm, 'Method', 'min-max', 'Range', [-1, 1]);
% Perform k-nearest neighbors imputation
%help knnimpute
%Ximp = knnimpute(Xnorm, 5, 'TreatMissing', 'ads', 'Distance', 'euc');
%
% CLUSTERING
addpath('../toolbox/clustering');
% Cluster data
%help kcentroids
[L, C, sumd] = kcentroids(Xnorm, 5, 'Treatmissing', 'ads', 'Distance', 'euc', ...
    'Replicates', 100, 'InitCrit', 'kmeans++', 'Start', []);
% Visualize clustering result
scatter_result(Xnorm, C, L);
%
% CLUSTER VALIDATION
addpath('../toolbox/cluster_validation');
% Use iterative clustering function to obtain centroids and labels for
% a selected range of K ([2, maxClusters])
[centers, labels] = iterative_kcentroids(Xm, 'MaxClusters', 10, 'TreatMissing', ... 
'exp', 'Distance', 'euc', 'Replicates', 100, 'InitCrit', 'kmeans++', 'Pipelined', false, ...
'UsePrevCent', true, 'ShowProgression', true);
% Define the used cluster validation indices
%help cluster_validation
indices = {@CalinskiHarabasz; @DaviesBouldin; @kCE; @PBM; @RayTuri; ... 
                @Silhouette; @WB; @WemmertGancarski};
% Optional way to define the indices
% indices = [{@CalinskiHarabasz, 'sqe'}; {@DaviesBouldin, 'euc'}; {@kCE, 'sqe'}; ... 
%             {@PBM, 'euc'}; {@RayTuri, 'sqe'}; {@Silhouette, 'euc'}; {@WB, 'sqe'}; ... 
%             {@WemmertGancarski, 'euc'}];
% Perform actual cluster validation
indices_values = cluster_validation(Xnorm, centers, labels, 'TreatMissing', 'exp', ... 
    'Distance', 'euc', 'ClusterInd', indices);
%
% Visualize the results
plot_indices(indices, indices_values);
%
% Restore default path
restoredefaultpath;

