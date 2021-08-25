# Nonrelational Network
 This repository contains a collection of code written in MATLAB. The goal of this code is to simulate spike trains using a Leaky Integrate-and-Fire model of firing with a nonrelational network structure of neurons with clustering.
 
 ## lif_network_postrotation.m
 This program contains code to initialize parameters, run network creation using create_clusters.m, run an LIF model using lif_sra_calculator_postrotation.m, and perform a number of post-run analyses.
 
 ## visualize_cluster_sequences.m
 This program contains code to visualize the sequence of clusters the spike sequences move through. The code generates not only single timestep visualizations but moving summed bin visualizations, and both can be seen in both regular weight and normalized form.
 
 ## calculate_plot_cluster_overlap.m
 This program calculates whether neurons firing at similar timepoints to each other are primarily from the same cluster, or different clusters. To do so, we use a sliding bin across the firing data and determine the average fraction of neurons in the same cluster(s) in each bin. If the bin has only a single neuron spiking, or no spikes, it will be left out of the calculation. If the bin has neurons from multiple different clusters, the calculation will be the sum of fraction of spiking neurons firing for each cluster represented, divided by the total number of clusters represented - getting an averge overlap calculation.
 
 ## Functions:
 
 ### lif_sra_calculator_postrotation.m
 This function uses the leaky integrate-and-fire model of  neuronal firing to calculate the trajectory of membrane potentials, currents, etc... that take place in a particular network with a particular set of parameters and initialization.
 
 ### create_clusters.m
 This function generates the network clusters and connections based on the number of neurons, number of clusters, number of neurons per cluster, and the probability of connections within a cluster.
 
 ### calculate_trajectory_similarity.m
 This function calculates a number of metrics of similarity between firing sequences. Specifically, it calculates Spearman's rank correlation rhos for sequences including and excluding nonfiring neurons from the ranks.
 
 ### calculate_cluster_overlap.m
 This function calculates the overlap of which clusters firing neurons belong to. For each event, using a sliding bin it calculates the fraction of overlap of all firing neurons for the full set of clusters they belong to.
