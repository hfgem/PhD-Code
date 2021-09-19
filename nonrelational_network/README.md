# Nonrelational Network
 This repository contains a collection of code written in MATLAB. The goal of this code is to simulate spike trains using a Leaky Integrate-and-Fire model of firing with a nonrelational network structure of neurons with clustering.
 
 ## lif_network_postrotation.m
 This program contains code to initialize parameters, run network creation using create_clusters.m, run an LIF model using lif_sra_calculator_postrotation.m, and perform a number of post-run analyses.
 
 ## viualize_cluster_sequences.m
 This program visualizes a sequences of clusters a spike sequences progresses through, based on outputs from lif_network_postrotation.m
 
 ## visualize_voltage_traces.m
 This program visualizes how a spike in one neuron affects connected neurons by looking at membrane potential deflections.

 ## Functions:
 
 ### calculate_trajectory_similarity_mi.m
 This function calculates the Matching Index (MI) (from Vas et al.) of trajectory similarity between firing sequences. Specifically, it calculates the MI for neuron firing orders including nonspiking neurons at the end of the order, as well as excluding nonspiking neurons.
 
 ### calculate_trajectory_similarity_spearmans.m
 This function calculates the Spearman's Rank Correlation Index of trajectory similarity between firing sequences with unique ranks. Specifically, it calculates the index for neuron firing orders including nonspiking neurons at the end of the order, as well as excluding nonspiking neurons.
 
 ### comp_percentile.m
 This function computes a percentile of a value against a dataset at a 10^(-2) accuracy.
 
 ### create_clusters.m
 This function generates the network clusters and connections based on the number of neurons, number of clusters, number of neurons per cluster, and the probability of connections within a cluster.
 
 ### generate_shuffled_trajectories.m
 This function generates shuffled firing sequences based on the statistics of real trajectories from network simulations.
 
 ### lif_sra_calculator_postrotation.m
 This function uses the leaky integrate-and-fire model of  neuronal firing to calculate the trajectory of membrane potentials, currents, etc... that take place in a particular network with a particular set of parameters and initialization.
 
 
 
 
 
 
