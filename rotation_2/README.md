# Rotation 2
 This repository contains code written in MATLAB during my second rotation. The goal is to analyze protein localization in Drosophila melanogaster third instar larval neuromuscular junctions (NMJs).
 
 ## NMJ_Segmentation_Analysis.m
 This program automatically loads and analyses data outputs from the PAZ Analysis ImageJ Macro developed by Dr. Steven Del Signore of the Rodal Lab at Brandeis University. The program takes inputs of which data measures are of interest and analyzes the data for whether it comes from a normal distribution and then tests if the control and experimental results come from the same or different distributions (these results are output in 'results.mat'). It also creates combined and separate box plots of the data measures desired.