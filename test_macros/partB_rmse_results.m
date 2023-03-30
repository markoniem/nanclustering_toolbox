% Description: 
% Evaluate performance of K-spatial-medians clustering based on Euclidean 
% distance estimation agains K-spatial-medians clustering based on available
% data strategy. Data sets consist of 20 % missing values. RMSEs values are 
% computed between real centroids and centroids obtained with clustering
% algorithms.
%
% Results are available in a xlsx file.  
%
% Inputs: 
%       See parameters function
%
% Output:
%       Computed RMSE errors between real centroids and clustered centroids
%
clear, clc, close all;
addpath('../datasets');
addpath('../toolbox/general');
addpath('../toolbox/clustering');
load('realCenters');
%
params = parameters();
datasets = params.datasets;
correctResults = params.correctResults;
probmiss = params.probmiss;
distance = params.distance;
replicates = params.replicates;
repetitions = params.repetitions;
initcrit = params.initcrit;
%
centroids = cell(length(datasets),repetitions,3);
labels = cell(length(datasets),repetitions,3);
realLabelsArr = cell(length(datasets),repetitions);
realCentersArr = cell(length(datasets),repetitions);
for i = 1:length(datasets)
    fprintf('Dataset: %s\n',datasets{i});
    for j = 1:repetitions
        fprintf('Repetition: %d\n',j);
        dataset = datasets{i};
        load(dataset);
        X = normalizedata(X,'Method','min-max','Range',[-1, 1]);
        [Xm, I] = genmissdata(X,probmiss);
        Xm = Xm(I,:);
        %
        % Update locations of real centers and save result 
        [Lreal, Creal, ~] = kcentroids(Xm,correctResults(i),'TreatMissing','ads', ...
            'Distance',distance,'Replicates',replicates,'InitCrit','kmeans++','Start',realCenters{i});       
        realCentersArr{i,j} = Creal;
        realLabelsArr{i,j} = Lreal;
        % K-spatial-medians clustering based on available data strategy
        [Lpred1, Cpred1, ~] = kcentroids(Xm,correctResults(i),'TreatMissing','ads', ...
            'Distance',distance,'Replicates',replicates,'InitCrit','kmeans++');
        [Cpred1, Lpred1] = correctOrdering(Creal, Cpred1, Lpred1);
        % K-spatial-medians clustering based on expected euclidean
        % distance strategy
        [Lpred2, Cpred2, ~] = kcentroids(Xm,correctResults(i),'TreatMissing','exp', ...
            'Distance',distance,'Replicates',replicates,'InitCrit','kmeans++');
        [Cpred2, Lpred2] = correctOrdering(Creal, Cpred2, Lpred2);        
        % Pipeline results for K-spatial-medians with available data
        % strategy
        [Lpred3, Cpred3, ~] = kcentroids(Xm,correctResults(i),'TreatMissing','ads', ...
            'Distance',distance,'Replicates',1,'InitCrit','kmeans++','Start',Cpred2);    
        [Cpred3, Lpred3] = correctOrdering(Creal, Cpred3, Lpred3);        
        centroids{i,j,1} = Cpred1;
        centroids{i,j,2} = Cpred2;
        centroids{i,j,3} = Cpred3;
        labels{i,j,1} = Lpred1;
        labels{i,j,2} = Lpred2;
        labels{i,j,3} = Lpred3;        
    end
end
save('centroids.mat','centroids');
save('realCentersArr.mat','realCentersArr');
%
%
rmse_errors = zeros(length(datasets),repetitions,4);
%
% Compute RMSE error between real centers and clustered real centers
% obtained with data sets consisting 20% missing values.
for i = 1:length(datasets)
    Creal = realCenters{i};
    alpha = 1/size(Creal,2);
    for j = 1:repetitions
        Cpred = realCentersArr{i,j};
        rmse = 0;   
        for k = 1:size(Creal,1)
            rmse = rmse + sqrt(alpha*sum((Creal(k,:)-Cpred(k,:)).^2));    
        end
        rmse = rmse / size(Creal,1);
        rmse_errors(i,j,1) = rmse; 
    end
end
%
% Compute RMSE error between real centers and clustering results based on
% available data distance computation strategy. Data sets consist of 
% 20% missing values.
for i = 1:length(datasets)
    Creal = realCenters{i};
    alpha = 1/size(Creal,2);
    for j = 1:repetitions
        Cpred = centroids{i,j,1};
        rmse = 0;   
        for k = 1:size(Creal,1)
            rmse = rmse + sqrt(alpha*sum((Creal(k,:)-Cpred(k,:)).^2));    
        end
        rmse = rmse / size(Creal,1);
        rmse_errors(i,j,2) = rmse; 
    end
end
%
% Compute RMSE error between real centers and clustering results based on
% expected euclidean distance computation strategy. Data sets consist of 
% 20% missing values.
for i = 1:length(datasets)
    Creal = realCenters{i};
    alpha = 1/size(Creal,2);
    for j = 1:repetitions
        Cpred = centroids{i,j,2};
        rmse = 0;   
        for k = 1:size(Creal,1)
            rmse = rmse + sqrt(alpha*sum((Creal(k,:)-Cpred(k,:)).^2));    
        end
        rmse = rmse / size(Creal,1);
        rmse_errors(i,j,3) = rmse; 
    end
end
%
% Compute RMSE error between real centers and clustering results based on
% pipelined expected euclidean - available data distance computation strategy. 
% Data sets consist of 20% missing values.
for i = 1:length(datasets)
    Creal = realCenters{i};
    alpha = 1/size(Creal,2);
    for j = 1:repetitions
        Cpred = centroids{i,j,3};
        rmse = 0;   
        for k = 1:size(Creal,1)
            rmse = rmse + sqrt(alpha*sum((Creal(k,:)-Cpred(k,:)).^2));    
        end
        rmse = rmse / size(Creal,1);
        rmse_errors(i,j,4) = rmse; 
    end
end
errors1 = mean(rmse_errors,2);
stds1 = std(rmse_errors,0,2);
%
rmse_errors = zeros(length(datasets),repetitions,3);
%
% Compute RMSE error between clustering results obtained using real centers 
% and randomly selected centers (utilizing kmeans++ algorithm) in 
% initialization of K-spatial-medians clustering that uses available data 
% distance computation. Data sets consits of 20% missing values.
for i = 1:length(datasets)
    for j = 1:repetitions
        Creal = realCentersArr{i,j};
        Cpred = centroids{i,j,1};
        alpha = 1/size(Creal,2);
        rmse = 0;   
        for k = 1:size(Creal,1)
            rmse = rmse + sqrt(alpha*sum((Creal(k,:)-Cpred(k,:)).^2));    
        end
        rmse = rmse / size(Creal,1);
        rmse_errors(i,j,1) = rmse; 
    end
end
%
% Compute RMSE error between clustering results obtained using real centers 
% and randomly selected centers (utilizing kmeans++ algorithm) in 
% initialization of K-spatial-medians clustering. Real centers are clustered 
% based on available data strategy distance computation and randomly selected 
% centers are clustered based on expected euclidean distance strategy. 
% Data sets consits of 20% missing values.
for i = 1:length(datasets)
    for j = 1:repetitions
        Creal = realCentersArr{i,j};
        Cpred = centroids{i,j,2};
        alpha = 1/size(Creal,2);
        rmse = 0;   
        for k = 1:size(Creal,1)
            rmse = rmse + sqrt(alpha*sum((Creal(k,:)-Cpred(k,:)).^2));    
        end
        rmse = rmse / size(Creal,1);
        rmse_errors(i,j,2) = rmse; 
    end
end
%
% Compute RMSE error between clustering results obtained using real centers 
% and randomly selected centers (utilizing kmeans++ algorithm) in 
% initialization of K-spatial-medians clustering. Real centers are clustered 
% based on available data strategy distance computation and randomly selected 
% centers are clustered based on pipelined expected euclidean - available data 
% clustering strategy. Data sets consits of 20% missing values.
for i = 1:length(datasets)
    for j = 1:repetitions
        Creal = realCentersArr{i,j};
        Cpred = centroids{i,j,3};
        alpha = 1/size(Creal,2);
        rmse = 0;   
        for k = 1:size(Creal,1)
            rmse = rmse + sqrt(alpha*sum((Creal(k,:)-Cpred(k,:)).^2));    
        end
        rmse = rmse / size(Creal,1);
        rmse_errors(i,j,3) = rmse; 
    end
end
errors2 = mean(rmse_errors,2);
stds2 = std(rmse_errors,0,2);
%
% Save results to .xlsx file
%
error_arr1 = cell(4,length(datasets)); 
error_arr2 = cell(3,length(datasets)); 
for i = 1:4
    for j = 1:length(datasets)
        error_arr1{i,j} = sprintf('%.3f (%.3f)',errors1(j,:,i),stds1(j,:,i));
    end
end
for i = 1:3
    for j = 1:length(datasets)
        error_arr2{i,j} = sprintf('%.3f (%.3f)',errors2(j,:,i),stds2(j,:,i));
    end
end
for i = 1:length(datasets), datasets{i} = sprintf('%s',datasets{i}); end
rows1 = cell(4,1);
rows1{1,1} = 'ADS: Real Cent. Init'; 
rows1{2,1} = 'ADS';
rows1{3,1} = 'EED';
rows1{4,1} = 'EED-ADS';
rows2 = cell(3,1); 
rows2{1,1} = 'ADS';
rows2{2,1} = 'EED';
rows2{3,1} = 'EED-ADS';
for i = 1:2
    if i == 1
        R = error_arr1;
        T = table(rows1,R);
        resultsTable = splitvars(T);
        resultsTable.Properties.VariableNames = [' ', datasets'];
        writetable(resultsTable,'partB_RMSE_Results.xlsx', 'Sheet', sprintf('Sheet %d',i));
    else
        R = error_arr2;
        T = table(rows2,R);
        resultsTable = splitvars(T);
        resultsTable.Properties.VariableNames = [' ', datasets'];
        writetable(resultsTable,'partB_RMSE_Results.xlsx', 'Sheet', sprintf('Sheet %d',i));        
    end
end
%
restoredefaultpath;


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
params.probmiss = 0.20;
params.distance = 'euc';
params.replicates = 100;
params.repetitions = 100;
params.initcrit = 'kmeans++';

end


