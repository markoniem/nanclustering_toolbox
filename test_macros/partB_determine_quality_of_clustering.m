% Description: 
% Uses external cluster validation indices to determine quality of clustering. 
% Saves validation indices results to xlsx file.
%
% Inputs: 
%       See parameters function
%
% Outputs:
%       indices_values - Values of external validation indice 
%
clear, clc, close all;
addpath('../datasets');
addpath('../toolbox/general');
addpath('../toolbox/clustering');
addpath('../toolbox/cluster_validation/external');
load('csai');
load('real_centers')
%
params = parameters();
datasets = params.datasets;
correctResults = params.correctResults;
probmiss = params.probmiss;
clusterInd = params.clusterInd;
%
indicesValues = zeros(length(datasets),length(clusterInd),length(probmiss),2);
for i = 1:length(datasets)
    fprintf('Data set: %s\n',datasets{i})
    for j = 1:length(clusterInd)
        for k = 1:length(probmiss)
            Xm = datamatrices{i,k};
            Creal = realCenters{i};
            % Compute real cluster labels
            [Lreal, Creal, ~] = kcentroids(Xm,correctResults(i),'TreatMissing','ads','Distance', ...
                'euc','Replicates',1,'InitCrit','kmeans++','Start',Creal);
            % Compute predicted cluster labels using two different clustering methods
            %
            [Lpred1, Cpred1, ~] = kcentroids(Xm,correctResults(i),'TreatMissing','ads','Distance', ...
                'euc','Replicates',100,'InitCrit','kmeans++','Start',[]);
            %
            [~, Cpred2, ~] = kcentroids(Xm,correctResults(i),'TreatMissing','exp','Distance', ...
                'euc','Replicates',100,'InitCrit','kmeans++','Start',[]);
            % Pipeline clustering results based on expected distances for clustering based on available data strategy
            [Lpred2, Cpred2, ~] = kcentroids(Xm,correctResults(i),'TreatMissing','ads','Distance', ...
                'euc','Replicates',1,'InitCrit','kmeans++','Start',Cpred2);
            %
            % Make sure correct ordering of cluster centroids
            [~, Lpred1] = correctOrdering(Creal, Cpred1, Lpred1);
            [~, Lpred2] = correctOrdering(Creal, Cpred2, Lpred2);
            %
            clusterIndex = clusterInd{j};
            indValues1 = clusterIndex(Lreal,Lpred1);
            indValues2 = clusterIndex(Lreal,Lpred2);
            indicesValues(i,j,k,1) = indValues1;
            indicesValues(i,j,k,2) = indValues2;
        end
    end
end
fprintf('Done!\n');
restoredefaultpath;
%
% Save results to .xlsx file
%
externalIndchar = cell(length(clusterInd),1);
externalIndchar2 = cell(length(clusterInd),1);
probMisschar = cell(length(probmiss),1);
for i = 1:length(clusterInd), externalIndchar{i} = sprintf('%s',char(clusterInd{i,1})); end
for i = 1:length(clusterInd), externalIndchar2{i} = sprintf('%s_2',char(clusterInd{i,1})); end
for i = 1:length(probmiss), probMisschar{i} = sprintf('%s%%',string(100*probmiss(i))); end
for i = 1:length(probmiss)
    R1 = indicesValues(:,:,i,1);
    R2 = indicesValues(:,:,i,2);
    R = [R1, R2];
    T = table(datasets,R);
    resultsTable = splitvars(T);
    resultsTable.Properties.VariableNames = [' ', [externalIndchar', externalIndchar2']];
    writetable(resultsTable,'partB_Results3.xlsx', 'Sheet', probMisschar{i});
end

function [C, L] = correctOrdering(CReal, CPred, LPred)
% Description:
% Orders cluster centroids to same order in different matrices.  
%
% Inputs:
%         CReal - First set of cluster centroids
%         CPred - Second set of cluster centroids 
%
% Output:
%         C - Predicted cluster centroids in same order as in CReal  
%         L - Labels corresponding clusters in C  
%
I = zeros(size(CReal,1),1);
C = CPred;
L = zeros(size(LPred));
for m = 1:size(CPred,1)
    D = pdist2(CPred,CReal);
    [~, index] = min(D(:));
    [row, col] = ind2sub(size(D), index);
    CReal(col,:) = Inf(size(CReal(col,:)));
    CPred(row,:) = Inf(size(CPred(row,:)));
    I(col) = row;
end
C = C(I,:);
for m = 1:size(C,1)
    L(LPred==I(m)) = m;
end

end

function params = parameters()
% Description: 
% Return initial parameters for test macro. 
%
% Output:
%     params - Parameters of test macro
%
params.datasets = {'S1'; 'S2'; 'S3'; 'S4'; 'S5D2'; 'S2D2'; 'O200'; 'O2000'};
params.correctResults = [15, 15, 15, 15, 5, 2, 5, 5];
params.probmiss = [0.00; 0.05; 0.10; 0.20];
params.clusterInd = {@ACC; @ARI; @NMI};

end

