# Entropy Clustering
This project addresses the problem of density-based clustering where clusters have variant densities and arbitrary shapes. Regarding the hardness of parameter regularization in existing methods, an entropy-based framework is given to generate clusters in two stages. Motivated by the fact that two clusters with ambiguous boundaries might be hardly discriminated, it presents the concept of density reachability to restrict the distance-based connectivity. By measuring entropy via k-NN distances in the first stage, a density order can be determined to glue the density-dominated cores into their adjacent dominators w.r.t. the density estimations. It highlights a heuristic algorithm in the second stage to effectively extract density peaks, such that the candidate peaks could enjoy a nice priority of low intra-cluster and high inter-cluster distances. The experimental study is conducted on dozens of data sets. The comparison results demonstrate that with simple tuning parameters, the proposed method shows substantial superiority to existing density-based clustering algorithms especially when the clusters exhibit various densities and arbitrary shapes.
## Invovled competitors
DBSCAN
DPC
KMEANS
## Results
Dataset|Metric| KMEANS| DBSCAN| PARMS |DPC| PARMS| EC| PARMS|
--- | --- | --- | --- |--- |---| --- | --- | ---  
    |ACC |0.629 |0.977 |      |0.902 |      |0.992 |
Jain|AMI |0.556 |0.973 |1.456 |0.878 |0.368 |0.992 |3
    |ARI |0.683 |0.943 |/3    |0.844 |      |0.981 |
