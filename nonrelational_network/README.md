# Nonrelational Network
 This repository contains a collection of code written in MATLAB. The goal of this code is to simulate spike trains using a Leaky Integrate-and-Fire model of firing with a nonrelational network structure of neurons with clustering.
 
 ## lif_network_postrotation.m
 This program contains code to initialize parameters, run network creation using create_clusters.m, run an LIF model using lif_sra_calculator_postrotation.m, and perform a number of post-run analyses.
 
 ## lif_sra_calculator_postrotation.m
 This function uses the leaky integrate-and-fire model of  neuronal firing to calculate the trajectory of membrane potentials, currents, etc... that take place in a particular network with a particular set of parameters and initialization.
 
 ## create_clusters.m
 This function generates the network clusters and connections based on the number of neurons, number of clusters, number of neurons per cluster, and the probability of connections within a cluster.
 
 ## calculate_trajectory_similarity.m
 This function calculates a number of metrics of similarity between firing sequences. Specifically, it calculates Spearman's rank correlation rhos for sequences including and excluding nonfiring neurons from the ranks.
