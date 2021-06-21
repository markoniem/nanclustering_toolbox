% Description: 
% Visualize the data transformation results
clear, clc, close all;
addpath('../../toolbox/preprocess');
addpath('../../benchmark_data/O-sets/');
addpath('../../benchmark_data/Sim-sets/');
addpath('../../benchmark_data/S-sets/');
load('S4');
figure;
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



