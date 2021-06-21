function params = params()
% Description: 
% Return initial parameters for test macros. See more detailed
% instructions of paramaters from original macros.  
%
% Output:
%     params - Parameters of test macro
%
params.datasets = {'S1'; 'S2'; 'S3'; 'S4'; 'S5D2'; 'S2D2'; 'O200'; 'O2000'};
params.probmiss = [0.00; 0.05; 0.10; 0.20];
params.maxClusters = 20;
params.clustMethod = @kcentroids_expected;
params.distance = 'euc';
params.replicates = 100;
params.useprevCent = true;
params.initcrit = 'kmeans++';
params.pipelined = true;
params.showProgression = true;
params.distFuncInd = 'ads';
params.clusterInd = {@CalinskiHarabasz; @DaviesBouldin; @DaviesBouldin2; ... 
   @GenDunn; @kCE; @PBM; @RayTuri; @Silhouette; @WB; @WemmertGancarski};
    
