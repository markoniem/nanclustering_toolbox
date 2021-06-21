% Description: 
% Visualize the data tranformation results 
clear, clc, close all;
addpath('../../benchmark_data/M_spheres');
load('K15_Nk400_M50_SPEHERES_dc0.9.mat');
addpath('M_Spheres2D');
addpath('../../toolbox/preprocess');
load('K15_Nk400_M50_SPEHERES_dc0.9.mat');
X = datamatrices{1,1};
figure;
hold on;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
axis square;
box on;
scatter(X(:,1),X(:,2),'.k');
[Xmapped, kp] = datasetmap(X);
figure;
hold on;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
axis square;
box on;
scatter(Xmapped(:,1),Xmapped(:,2),'.k');
scatter(Xmapped(kp,1),Xmapped(kp,2),'ro','filled');
%
% Color figures
% figure;
% hold on;
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% axis square;
% box on;
% kNumber = size(C,1);
% for k = 1:kNumber
%     scatter(X(labels==k,1),X(labels==k,2),'.','CData',[rand rand rand]);
% end
% scatter(X(kp,1),X(kp,2),'ro','filled');



