% Description: 
%
% Performs pairwise distance computation between two data sets. 
% Pairwise distances can be computed using available data, partial
% distance, expected euclidean distance, or expected squared 
% euclidean distance strategy.
%
clear, clc, close all;
rng(1);
addpath('../datasets');
addpath('../toolbox/general');
load('real_centers');
%
% Centers of data sets are ordered as S1, S2, S3, S4, Sim5D2, Sim2D2, O200, and O2000.
% Use centers of S4 data set.
C = realCenters{4};
%
% Load S4 data set.
load('S4');
%
% Normalize data to a range of [-1, 1] and generate missing values.
X = normalizedata(X,'Method','min-max','Range',[-1,1]);
[X, I] = genmissdata(X,0.10);
X = X(I,:);
%
% Compute distances between observations and centroids.
% help nanpdist2
D = nanpdist2(X,C,'Distance','euc','TreatMissing','esd');
%
% Identify closests centroids
[~, L] = min(D,[],2);
%
% Visualize result
figure;
hold on;
for i = 1:size(C,1)
    L1 = (L==i);
    scatter(X(L1,1),X(L1,2));
end
scatter(C(:,1),C(:,2),'ro','filled');


