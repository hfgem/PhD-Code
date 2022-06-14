# ICA
 This repository contains a collection of code written to perform independent component analysis (ICA) using kurtosis maximization. This is a very basic implementation that takes aspects of ICA from the following papers:
 - "Independent Component Analysis, a new concept?" - Pierre Comon 1994
 - "Independent Component Analysis: Algorithms and Applications" - Aapo Hyvärinen and Erkki Oja, 2000
 
 # ICA_examples.m
 This file generates mixtures from signals and tests ICA results using both the manually written functions and the MATLAB rica function.
 
 # manual_kurtosis_ICA_rowwise.m
 This function performs kurtosis ICA using the Gram-Schmidt Orthogonalization, estimating one independent component at a time.
 
 # manual_kurtosis_ICA.m
 This function performs kurtosis ICA using Symmetric Orthogonalization, estimating all independent components in parallel.