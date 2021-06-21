% Description: 
% Generates missing values to complete data set(s)
%
% Inputs (see current values from 'parameters.m' file):
%        datasets - Input data set(s)
%        probmiss - Probabilities of missing values generated to data set(s)
%
% Output:
%    datamatrices - Data sets which consist of predefined numbers of missing
%                   values. First cell is for different data sets and second 
%                   cell is for different percentages of missing values 
%                   generated to data.
%
clear, clc, close all;
rng(1); % For reproduction purpose
addpath('../../benchmark_data/O-sets/');
addpath('../../benchmark_data/Sim-sets/');
addpath('../../benchmark_data/S-sets/');
addpath('../../toolbox/preprocess');
params = params();
datasets = params.datasets;
probmiss = params.probmiss;
datamatrices = cell(length(datasets),length(probmiss));
fprintf('Generating missing values...\n');
for i = 1:length(datasets)
    for j = 1:length(probmiss)
        load(char(datasets(i)));
        X = normalizedata(X,'min-max');
        Xm = genmissdata(X,probmiss(j));
        datamatrices{i,j} = Xm;
    end
end
fprintf('Done!\n');
save('datamatrices','datamatrices');



