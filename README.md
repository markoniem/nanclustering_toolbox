#### TOOLBOX ARTICLE ####
[M. Niemelä, S. Äyrämö, and T. Kärkkäinen, "Toolbox for Distance Estimation and Cluster Validation 
on Data With Missing Values," *in IEEE Access, vol. 10,* pp. 352-367, 2022](https://ieeexplore.ieee.org/document/9656159)

#### ABSTRACT ####
*"Missing data are unavoidable in the real-world application of unsupervised machine learning,
and their nonoptimal processing may decrease the quality of data-driven models. Imputation is a common
remedy for missing values, but directly estimating expected distances have also emerged. Because
treatment of missing values is rarely considered in clustering related tasks and distance metrics have
a central role both in clustering and cluster validation, we developed a new toolbox that provides a
wide range of algorithms for data preprocessing, distance estimation, clustering, and cluster validation
in the presence of missing values. All these are core elements in any comprehensive cluster analysis
methodology. We describe the methodological background of the implemented algorithms and present
multiple illustrations of their use. The experiments include validating distance estimation methods against
selected reference methods and demonstrating the performance of internal cluster validation indices. The
experimental results demonstrate the general usability of the toolbox for the straightforward realization of
alternate data processing pipelines."*

#### SUPPORTED FUNCTIONALITIES ####
**Distance estimation methods in case of incomplete data:**  
Available data strategy (ADS), Partial distance strategy (PDS), Expected squared Euclidean distance (ESD), Expected Euclidean distance (EED)

**Pre-processing methods:**  
Min-max and z-score scalings, missing values generation, CC-/ICC-/k-nearest neighbors imputation (support for ADS, PDS, ESD, and EED), expectation maximization with maximum likelihood convergence

**Clustering algorithms:**  
K-means (support for ADS, PDS, ESD, and EED), K-medians (support for ADS and PDS), K-spatial-medians (support for ADS, PDS, ESD, and EED). **NOTE:** <ins>Support for weights added on 28 March, 2023</ins> 

**Internal cluster validation indices:**  
Calinski-Harabasz, Davies-Bouldin, Davies-Bouldin*, Generalized Dunn, kCE-index, Pakhira-Bandyopadhyay-Maulik, Ray-Turi, Silhouette, WB-index, Wemmert-Gancarski (support for ADS, PDS, ESD, and EED). **NOTE:** <ins>Support for weights added on 28 March, 2023</ins> 
  
**External cluster validation indices:**  
Accuracy-index, Adjusted Rand index, Normalized mutual information

**Parallel and GPU computations:**  
Clustering and cluster validation methods support parallelled and GPU computations

#### REFERENCES #### 
[See toolbox article](https://ieeexplore.ieee.org/document/9656159)

**Weighted clustering:**    
M. Saarela and T. Kärkkäinen, "Do country stereotypes exist in PISA? a clustering approach for large, sparse, and weighted data," *In Proceedings of the 8th International Conference on Educational Data Mining, International Educational Data Mining Society,* EDM, pages 156-163, 2015
 
