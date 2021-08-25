%Calculate and Plot Cluster Overlap

%This code calculates whether neurons firing at similar timepoints to
%each other are primarily from the same cluster, or different clusters.
%To do so, we use a sliding bin across the firing data and determine the
%average fraction of neurons in the same cluster(s) in each bin. If the bin
%has only a single neuron spiking, or no spikes, it will be left out of the
%calculation. If the bin has neurons from multiple different clusters, the
%calculation will be the sum of fraction of spiking neurons firing for each
%cluster represented, divided by the total number of clusters represented -
%getting an averge overlap calculation.

%Load data and clear unnecessary variables from memory
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
slashes = find(net_save_path == '/');
save_path = net_save_path(1:slashes(end));

load(strcat(net_save_path,'/network.mat'))
cluster_mat = network.cluster_mat; %matrix of neurons in each cluster
clear network %get rid of the large memory suck

load(strcat(save_path,'/parameters.mat'))
n = parameters.n;
V_th = parameters.V_th; %threshold membrane potential (V) - update this from the parameter txt file
clusters = parameters.clusters;
clear parameters %get rid of the large memory suck

load(strcat(net_save_path,'/network_spike_sequences.mat'))
events = {network_spike_sequences.events};
clear network_spike_sequences %get rid of the large memory suck

load(strcat(net_save_path,'/network_var.mat'))
[~, init] = size(events);
V_m_struct = struct;
for i = 1:init
    V_m = network_var(i).V_m;
    V_m_struct(i).V_m = V_m;
    V_m_struct(i).spikes_V_m = V_m >= V_th;
end    
clear i V_m network_var %get rid of the large memory suck

%Calculate cluster overlaps
bin_frac = 0.1; %percent of total time represented to bin together in cluster overlap calc
[cluster_prog] = calculate_cluster_overlap(init,V_m_struct,cluster_mat, ...
    bin_frac, events);
save(strcat(net_save_path,'/cluster_prog.mat'),'cluster_prog')

%Visualize overlap over time in regular-time data
colors = {'r','g','b','m','k'};
f = figure;
axes = [];
subplot_n = sqrt(init);
if floor(subplot_n)==subplot_n %assume non-infinite
    subplot_x = subplot_n;
    subplot_y = subplot_n;
else
    subplot_x = floor(subplot_n);
    subplot_y = ceil(subplot_n);
end
clear subplot_n
for i = 1:init
    [~,event_num] = size(cluster_prog(i).event);
    ax = subplot(subplot_x,subplot_y,i);
    axes(end+1) = ax; %#ok<SAGROW>
    hold on
    lg = {};
    for j = 1:event_num
        plot(cluster_prog(i).event(j).cluster_overlap,colors{j})
        lg(end+1) = {strcat('Event #',string(j))};
    end
    ylabel('Amount of Overlap')
    xlabel('Event Time Step')
    title(strcat('Initialization #',string(i)))
    legend(lg)
end
sgtitle('Cluster Overlap - Over Time')
linkaxes(axes)

%Create image save path
cluster_save_path = strcat(net_save_path,'/cluster_plots/');
if ~isfolder(cluster_save_path)
    mkdir(cluster_save_path);
end

savefig(f,strcat(cluster_save_path,'/cluster_overlap.fig'))
saveas(f,strcat(cluster_save_path,'/cluster_overlap.jpg'))