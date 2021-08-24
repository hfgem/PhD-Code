%LIF Network Rotation Project Expanded

%% Save Path

save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder'); %Have user input where they'd like the output stored

%% Initialization

%Neuron and cluster counts
n = 200; %number of neurons
clusters = round(n/20); %number of clusters of neurons (for small n round(n/5), for large n round(n/20)) 
cluster_n = round(n/5); %number of neurons in a cluster (for small n round(n/3), for large n round(n/5)) 

%Interaction constants
t_max = 2; %maximum amount of time (s)
dt = 0.1*10^(-3); %timestep (s)
t_steps = t_max/dt; %number of timesteps in simulation
tau_syn_E = 20*10^(-3); %AMPA/NMDA synaptic decay time constant (s)
tau_syn_I = 5*10^(-3); %GABA synaptic decay time constant (s)
tau_sra = 250*10^(-3); %spike rate adaptation time constant (s)
E_K = -80*10^(-3); %potassium reversal potential (V) %-75 or -80 mV
E_L = -70*10^(-3); %leak reversal potential (V) %-60 - -70 mV range
G_L = 25*10^(-9); %leak conductance (S) %10 - 30 nS range
C_m = 0.5*10^(-9); %total membrane capacitance (F) %Huge range from 0.1 - 100 pF
V_th = -50*10^(-3); %threshold membrane potential (V)
V_reset = -70*10^(-3); %reset membrane potential (V)
V_syn_E = 0; %synaptic reversal potential (excitatory)
V_syn_I = -70*10^(-3); %synaptic reversal potential (inhibitory) %generally -70 pr -80 mV
E_syn_E = V_syn_E*ones(n,1); %vector of the synaptic reversal potential for excitatory connections
E_syn_I = V_syn_I*ones(n,1); %vector of the synaptic reversal potential for inhibitory connections
%______Split del_G_syn______
del_G_syn_E = 10*10^(-9); %synaptic conductance step following spike (S)
del_G_syn_I = 10*10^(-9); %synaptic conductance step following spike (S)
%___________________________
del_G_sra = 200*10^(-9); %spike rate adaptation conductance step following spike %ranges from 1-200 *10^(-9) (S)
connectivity_gain = 1; %amount to increase or decrease connectivity by with each spike (more at the range of 1.002-1.005)
IEI = 0.5; %inter-event-interval (s) the elapsed time between spikes to count separate events
IES = ceil(IEI/dt); %inter-event-steps = the number of steps to elapse between spikes

%How spikes are initiated:
%'cluster' sets a cluster to threshold;
%'current' means spikes depend on an input current of I_in; 
%'neuron' sets a random selection of 2% of neurons to threshold
type = 'cluster'; %'current'; %'neuron';

%How many tests of different initializations to run
if strcmp(type,'cluster')
    test_val_max = clusters; %Ensures every cluster is initialized
else
    test_val_max = 10; %This value can be modified as you see fit
end   

%Adding an input current to all cells (one of the options must be uncommented)
x_in = [0:dt:t_max];
% %Rhythmic current input: (uncomment if desired)
% I_coeff = 0; %5.1*10^(-10); %set to 0 for no input current
% I_Hz = 1; %frequency of input - 1 Hz or 60 bpm, 4-8 Hz for theta rhythm
% I_in = I_coeff*(0.5 + 0.5*sin(I_Hz*2*pi*x_in)); %Approximately theta wave input current
%Noisy input: (uncomment if desired)
I_coeff = 0; %set to 0 for no input current
I_scale = 1*10^(-9); %sets the scale of the current input
I_in = I_coeff*randn(n,t_steps+1)*I_scale; %Generally use -0.5-0.5 nA stimulus

%Calculate connection probabilites
conn_prob = 0.08; %set a total desired connection probability
npairs = n*(n-1); %total number of possible neuron connections
nclusterpairs = cluster_n*(cluster_n - 1)*clusters; %total number of possible intra-cluster connections
cluster_prob = min(conn_prob*npairs/nclusterpairs,1); %intra-cluster connection probability
p_E = 0.75; %probability of an excitatory neuron
p_I = 1 - p_E; %probability of an inhibitory neuron
n_I = round(p_I*n); %number of inhibitory neurons
all_indices = [1:n];

%OPTIONAL: Save parameters to a .txt file (note, must have empty Workspace
%befor running this section for this to work)
param_names = who;
fileID = fopen(strcat(save_path,'/param_file.txt'), 'wt'); % Open and create file in text writing mode.
fprintf(fileID, 'Parameters:\r\n');
fprintf(fileID, ' =======================================\r\n\r\n\r\n');
for i = 1:length(param_names)
    value = eval(param_names{i});
    [a,b] = size(value);
    if strcmp(class(value),'double') && a*b == 1 %#ok<STISA> %Only saves doubles of single value
        fprintf(fileID, strcat(param_names{i},'\t%s\r\n'), sprintf('%0.5e',value)); %write the variable
    elseif strcmp(param_names{i},'type')
        fprintf(fileID, strcat(param_names{i},'\t%s\r\n'), value); %write the variable
    end
    clear value a b
end
clear i
fclose(fileID);

%% Create networks and test spike progressions
%Runs through a series of different random number generator seeds to change
%the network connectivity and setup, and then automatically outputs
%sequence data to a folder

for i = 1:10 %how many different network structures to test
    rng(i) %set random number generator for network structure
    
    %CREATE NETWORK SAVE PATH
    net_save_path = strcat(save_path,'/network_',string(i));
    if ~isfolder(net_save_path)
        mkdir(net_save_path);
    end
    
    %SET UP NETWORK
    [cluster_mat, conns] = create_clusters(n, clusters, cluster_n, ...
        cluster_prob);
    conns_copy = conns; %just a copy of the connections to maintain for reset runs if there's "plasticity"
    %Randomize excitatory and inhibitory connection strengths based on selected
    %probability.
    I_indices = randi(n,[n_I,1]); %indices of inhibitory neurons
    all_indices(I_indices) = 0;
    E_indices = find(all_indices)'; %indices of excitatory neurons
    n_EE = sum(conns(E_indices,E_indices),'all'); %number of E-E connections
    n_EI = sum(conns(E_indices,I_indices),'all'); %number of E-I connections
    n_II = sum(conns(I_indices,I_indices),'all'); %number of I-I connections
    n_IE = sum(conns(I_indices,E_indices),'all'); %number of I-E connections
    
    %SAVE NETWORK STRUCTURE
    network = struct;
    network(1).cluster_mat = cluster_mat;
    network(1).conns = conns;
    network(1).I_indices = I_indices;
    network(1).E_indices = E_indices;
    save(strcat(net_save_path,'/network.mat'),'network')
    clear network %to save space
    
    %RUN MODEL AND CALCULATE
    %Run through every cluster initialization and store relevant data and
    %calculations
    network_var = struct;
    network_spike_sequences = zeros(clusters,n);
    network_spike_ranks = zeros(clusters,n);
    non_spiking_neurons = zeros(clusters,n);
    network_cluster_sequences = struct; 
    
    for j = 1:test_val_max
        seed = j;
        
        %Create Storage Variables
        I_syn = zeros(n,t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
        %synaptic conductances for each neuron at each timestep
        G_syn_I = zeros(n,t_steps+1); %conductance for presynaptic inhibitory (S)
        G_syn_E = zeros(n,t_steps+1); %conductance for presynaptic excitatory (S)
        V_m = zeros(n,t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = V_reset + randn([n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
        G_sra = zeros(n,t_steps+1); %refractory conductance for each neuron at each timestep
        
        %Run model
        [V_m, G_sra, G_syn_I, G_syn_E, I_syn] = lif_sra_calculator_postrotation(n, ...
            cluster_mat,conns,V_m,V_reset,V_th,G_sra,del_G_sra,...
            G_syn_I, G_syn_E, del_G_syn_E,del_G_syn_I,I_syn,E_syn_E,E_syn_I,E_K,E_L,...
            G_L,C_m,dt,tau_sra,tau_syn_E,tau_syn_I,connectivity_gain,...
            I_indices,E_indices,I_in,seed,t_steps,type);
        network_var(j).V_m = V_m;
        network_var(j).G_sra = G_sra;
        network_var(j).G_syn_I = G_syn_I;
        network_var(j).G_syn_E = G_syn_E;
        network_var(j).I_syn = I_syn;
        
        %Find spike sequence
        spikes_V_m = V_m >= V_th;
        [spikes_x,spikes_t] = find(spikes_V_m);
        max_time = max(spikes_t);
        spike_order = unique(spikes_x,'stable');
        for k = 1:length(spike_order)
            n_ind = spike_order(k);
            network_spike_ranks(j,n_ind) = k;
        end    
        for k = 1:n %make sure all indices are accounted for (non-spiking neurons)
            if isempty(find(spike_order == k,1))
                spike_order(end+1) = k; %#ok<SAGROW>
                non_spiking_neurons(j,k) = 1;
            end
        end
        network_spike_sequences(j,:) = spike_order;
        clear k spikes_x
        
        %Find maximum firing rate + average maximum firing rates of neurons
        all_fr = sum(spikes_V_m,2)/t_max;
        max_fr = max(all_fr);
        avg_fr = mean(all_fr);
        %USE THIS AS AN INDICATOR OF GOOD SEQUENCES, RATHER THAN VISUALIZING
        %In hippocampus, Jadhav lab uses 3 Hz as a place-cell cutoff
        if and(max_fr <= 3, max_fr >= 1) || and(avg_fr <= 3, avg_fr >= 0.5)
            %Visualize re-ordered spike sequence
            reordered_spikes = zeros(size(spikes_V_m));
            reordered_V_m = zeros(size(V_m));
            for i = 1:n
                reordered_spikes(i,:) = spikes_V_m(spike_order(i),:);
                reordered_V_m(i,:) = V_m(spike_order(i),:);
            end
            f = figure;
            imagesc(reordered_spikes(:,1:max(spikes_t)))
            xticks(round(linspace(1,max(spikes_t),10)))
            xt = get(gca, 'XTick');
            xtlbl = round(linspace(0, t_max, numel(xt)),2);
            colormap(flip(gray))
            if strcmp(type,'cluster')
                title(strcat('Spiking Behavior: cluster = ',string(j)),'FontSize',16)
            else
                title('Spiking Behavior','FontSize',16)
            end
            xlabel('Time (s)','FontSize',16)
            ylabel('Reordered Neuron Number','FontSize',16)
            set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
            savefig(f,strcat(net_save_path,'/',type,'_',string(j),'firing_sequence.fig'))
            saveas(f,strcat(net_save_path,'/',type,'_',string(j),'firing_sequence.jpg'))
            close(f)
            clear i f

            %Find event times
            events = []; 
            last_start = spikes_t(1);
            last_time = spikes_t(1);
            for t_i = 2:length(spikes_t)
                s_i = spikes_t(t_i);
                if s_i - last_time <= IES
                    last_time = s_i; %update the last time to be this one
                else
                    events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last range of spikes to the events vector
                    last_start = s_i;
                    last_time = s_i;
                end
            end
            if last_start ~= last_time %if the last event is a single spike, we don't care about it
                events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last interval
            end
            [num_events,~] = size(events);
            
            %Find cluster sequence per event by moving bin
            bin_width = 5*10^(-3); %5 ms bin
            bin_size = ceil(bin_width/dt); %number of timesteps to use in a bin
            for e_i = 1:num_events
                cluster_spikes = cluster_mat*spikes_V_m(:,events(e_i,1):events(e_i,2));
                cluster_mov_sum = movsum(cluster_spikes',bin_size)';
                normalized_cluster_spikes = cluster_spikes ./ sum(cluster_spikes,1);
                normalized_cluster_spikes(isnan(normalized_cluster_spikes)) = 0;
                normalized_cluster_mov_sum = cluster_mov_sum ./ sum(cluster_mov_sum,1);
                normalized_cluster_mov_sum(isnan(normalized_cluster_mov_sum)) = 0;
                network_cluster_sequences(j).clusters.(strcat('cluster_',string(e_i))) = cluster_spikes;
                network_cluster_sequences(j).movsum.(strcat('cluster_',string(e_i))) = cluster_mov_sum;
                network_cluster_sequences(j).normalized_clusters.(strcat('cluster_',string(e_i))) = normalized_cluster_spikes;
            end
            clear bin_width bin_size cluster_spikes cluster_mov_sum e_i
        end
    end
    
    %SAVE NETWORK DATA
    save(strcat(net_save_path,'/network_var.mat'),'network_var','-v7.3')
    save(strcat(net_save_path,'/network_spike_sequences.mat'),'network_spike_sequences')
    save(strcat(net_save_path,'/network_spike_ranks.mat'),'network_spike_ranks')
    save(strcat(net_save_path,'/non_spiking_neurons.mat'),'non_spiking_neurons')
    save(strcat(net_save_path,'/network_cluster_sequences.mat'),'network_cluster_sequences')
end

%% Calculate trajectory similarity
%This code block does not require any of the earlier code blocks to run. It
%simply requires stored data from prior runs of earlier code blocks.

% %Select and load data to analyze
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
load(strcat(net_save_path,'/network_spike_sequences.mat'))
load(strcat(net_save_path,'/non_spiking_neurons.mat'))
load(strcat(net_save_path,'/network_spike_ranks.mat')) %ranks not including nonspiking

[inits,n] = size(network_spike_sequences); %Grab data sizes

[neuron_ranks, rho_unique, ranks, rho_nonunique, ranks_mod] = calculate_trajectory_similarity(n, ...
    inits, network_spike_sequences, network_spike_ranks, non_spiking_neurons);


average_rank = mean(ranks(ranks~=1),'all');
std_rank = std(ranks(ranks~=1),[],'all');
average_rank_mod = mean(ranks_mod(ranks_mod~=1),'all');
std_rank_mod = std(ranks_mod(ranks_mod~=1),[],'all');

%% Visualize cluster sequences
%This code block visualizes the sequence of clusters a particular spike
%sequence progresses through.

%Select and load data to analyze
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
load(strcat(net_save_path,'/network_spike_sequences.mat'))
load(strcat(net_save_path,'/network_cluster_sequences.mat'))
load(strcat(net_save_path,'/network_var.mat'))

V_m = network_var(1).V_m;
spikes_V_m = V_m >= -50*10^(-3);

%Grab relevant information
[clusters,~] = size(network_spike_sequences); %Grab data sizes

%Visualize cluster sequences
f3 = figure;
subplot_n = sqrt(clusters);
if floor(subplot_n)==subplot_n %assume non-infinite
    subplot_x = subplot_n;
    subplot_y = subplot_n;
else
    subplot_x = floor(subplot_n);
    subplot_y = ceil(subplot_n);
end
clear subplot_n
%Movsum Cluster Plots
axes = [];
for k = 1:clusters
    ax = subplot(subplot_y,subplot_x,k);
    axes(end+1) = ax; %#ok<SAGROW>
    imagesc(network_cluster_sequences(k).movsum)
    title(strcat('Moving Sum Cluster Sequence for #',string(k)))
    xlabel('Binned Spike Times')
    ylabel('Cluster Number')
end
linkaxes(axes)
savefig(f3,strcat(net_save_path,'/movsum_cluster_sequence.fig'))
%close(f2)
clear k ax axes tick_vals
%Regular Cluster Plots
f3 = figure;
axes = [];
for k = 1:clusters
    ax = subplot(subplot_y,subplot_x,k);
    axes(end+1) = ax; %#ok<SAGROW>
    imagesc(network_cluster_sequences(k).clusters)
    title(strcat('Cluster sequence for #',string(k)))
    xlabel('Spike Times')
    ylabel('Cluster Number')
end
linkaxes(axes)
savefig(f3,strcat(net_save_path,'/cluster_sequence.fig'))
%close(f3)
clear k ax axes tick_vals

%% Calculate Cluster Firing Correlation
%This code block calculates whether neurons firing at similar timepoints to
%each other are primarily from the same cluster, or different clusters.
%To do so, we use a sliding bin across the firing data and determine the
%average fraction of neurons in the same cluster(s) in each bin. If the bin
%has only a single neuron spiking, or no spikes, it will be left out of the
%calculation. If the bin has neurons from multiple different clusters, the
%calculation will be the sum of fraction of spiking neurons firing for each
%cluster represented, divided by the total number of clusters represented -
%getting an averge overlap calculation.
save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Parameter Save Folder'); %Have user select folder where parameter file is saved
net_save_path = uigetdir(save_path,'Select Network Save Folder'); %Have user input the specific network
load(strcat(net_save_path,'/network_var.mat'))
load(strcat(net_save_path,'/network.mat'))

V_th = -50*10^(-3); %threshold membrane potential (V) - update this from the parameter txt file

%Grab relevant information and clear excess
[~, init] = size(network_var);
V_m_struct = struct;
for i = 1:init
    V_m = network_var(i).V_m;
    V_m_struct(i).V_m = V_m;
    V_m_struct(i).spikes_V_m = V_m >= V_th;
end    
clear i V_m network_var %get rid of the large memory suck
[n, t] = size(V_m_struct(1).V_m); %number of neurons x number of timesteps
cluster_mat = network.cluster_mat; %matrix of neurons in each cluster
[clusters, ~] = size(cluster_mat); %number of clusters
clear network %get rid of the large memory suck

%Calculate cluster overlaps
bin_frac = 0.005; %percent of total time represented to bin together in cluster overlap calc
cluster_prog = struct; %create structure to save all data
for i = 1:init
    %Pull spike data
    spikes_V_m = V_m_struct(i).spikes_V_m;
    [~,spikes_t] = find(spikes_V_m);
    total_spikes = length(unique(spikes_t));
    cluster_prog(i).total_spikes = total_spikes;
    if total_spikes > 1 %First make sure we have multiple spikes to work with
        collapsed_spikes_V_m = spikes_V_m(:,unique(spikes_t)); %Only look at spike times

        %First find bin size for each dataset
        bin_size = ceil(bin_frac*t); %how many timesteps to bin together
        collapsed_bin_size = 1; %ceil(bin_frac*total_spikes); %how many collapsed timesteps to bin together

        %Run through each dataset and calculate overlap fractions
        cluster_number = [];
        neuron_number = [];
        cluster_overlap = [];
        collapsed_cluster_number = [];
        collapsed_neuron_number = [];
        collapsed_cluster_overlap = [];
        for j = 1:t-bin_size %regular data
            spike_cut = spikes_V_m(:,j:j+bin_size); %n x bin_size
            n_spiking = sum(sum(spike_cut,2) > 0); %number of neurons that fire
            c_spiking = cluster_mat*(sum(spike_cut,2) > 0); %number of neurons per cluster that fire
            n_c = sum(c_spiking > 0); %number of clusters
            if n_spiking > 1 %only look at bins with more than 1 spike
                if n_spiking*n_c > 0 %leave out NaN bins
                    avg_overlap = sum(c_spiking)/(n_spiking*n_c);
                    cluster_overlap = [cluster_overlap, avg_overlap]; %#ok<AGROW>
                    cluster_number = [cluster_number, n_c]; %#ok<AGROW>
                    neuron_number = [neuron_number, n_spiking]; %#ok<AGROW>
                end
            end
        end
        clear j spike_cut n_spiking c_spiking n_c avg_overlap
        for j = 1:total_spikes - 1 - collapsed_bin_size %collapsed data
            spike_cut = collapsed_spikes_V_m(:,j:j+collapsed_bin_size);
            n_spiking = sum(sum(spike_cut,2) > 0); %number of neurons that fire
            c_spiking = cluster_mat*(sum(spike_cut,2) > 0); %gives the number of neurons per cluster that fire
            n_c = sum(c_spiking > 0); %number of clusters
            if n_spiking > 1 %only look at bins with more than 1 spike
                if n_spiking*n_c > 0 %leave out NaN bins
                    avg_overlap = sum(c_spiking)/(n_spiking*n_c);
                    collapsed_cluster_overlap = [collapsed_cluster_overlap, avg_overlap]; %#ok<AGROW>
                    collapsed_cluster_number = [collapsed_cluster_number, n_c];  %#ok<AGROW>
                    collapsed_neuron_number = [collapsed_neuron_number, n_spiking]; %#ok<AGROW>
                end
            end
        end
        clear j spike_cut n_spiking c_spiking n_c avg_overlap
        
        %Store relevant information
        cluster_prog(i).bin_size = bin_size;
        cluster_prog(i).collapsed_bin_size = collapsed_bin_size;
        cluster_prog(i).cluster_number = cluster_number;
        cluster_prog(i).avg_cluster_number = mean(cluster_number);
        cluster_prog(i).neuron_number = neuron_number;
        cluster_prog(i).avg_neuron_number = mean(neuron_number);
        cluster_prog(i).cluster_overlap = cluster_overlap;
        cluster_prog(i).avg_cluster_overlap = mean(cluster_overlap);
        cluster_prog(i).collapsed_cluster_number = collapsed_cluster_number;
        cluster_prog(i).avg_collapsed_cluster_number = mean(collapsed_cluster_number);
        cluster_prog(i).collapsed_neuron_number = collapsed_neuron_number;
        cluster_prog(i).avg_collapsed_neuron_number = mean(collapsed_neuron_number);
        cluster_prog(i).collapsed_cluster_overlap = collapsed_cluster_overlap;
        cluster_prog(i).avg_collapsed_cluster_overlap = mean(collapsed_cluster_overlap);
    end
    
    %Clear variables for space
    clear spikes_V_m spikes_x spikes_t collapsed_spikes_V_m bin_size ...
        collapsed_bin_size
end
clear i
save(strcat(net_save_path,'/cluster_prog.mat'),'cluster_prog')

%Visualize overlap over time in regular-time data
subplot_n = sqrt(init);
if floor(subplot_n)==subplot_n %assume non-infinite
    subplot_x = subplot_n;
    subplot_y = subplot_n;
else
    subplot_x = floor(subplot_n);
    subplot_y = ceil(subplot_n);
end
clear subplot_n
axes = [];
f = figure;
hold on
for i = 1:init
    ax = subplot(subplot_y,subplot_x,i);
    axes(end+1) = ax; %#ok<SAGROW>
    plot(cluster_prog(i).cluster_overlap,'b')
    yline(cluster_prog(i).avg_cluster_overlap,'r')
    legend('Trend','Average')
    title(strcat('Initialization #',string(i)))
end   
sgtitle('Cluster Overlap - Over Time')
linkaxes(axes)
savefig(f,strcat(net_save_path,'/cluster_overlap.fig'))
saveas(f,strcat(net_save_path,'/cluster_overlap.jpg'))

%Visualize overlap over time in collapsed data
subplot_n = sqrt(init);
if floor(subplot_n)==subplot_n %assume non-infinite
    subplot_x = subplot_n;
    subplot_y = subplot_n;
else
    subplot_x = floor(subplot_n);
    subplot_y = ceil(subplot_n);
end
clear subplot_n
axes = [];
f2 = figure;
hold on
for i = 1:init
    ax = subplot(subplot_y,subplot_x,i);
    axes(end+1) = ax; %#ok<SAGROW>
    plot(cluster_prog(i).collapsed_cluster_overlap,'b')
    yline(cluster_prog(i).avg_collapsed_cluster_overlap,'r')
    legend('Trend','Average')
    title(strcat('Initialization #',string(i)))
end   
sgtitle('Collapsed Cluster Overlap - Over Time')
linkaxes(axes)
savefig(f2,strcat(net_save_path,'/collapsed_cluster_overlap.fig'))
saveas(f2,strcat(net_save_path,'/collapsed_cluster_overlap.jpg'))


%% Visualize network structure
%Select data to visualize
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
load(strcat(net_save_path,'/network.mat'))

cluster_mat = network.cluster_mat;

%First create indices for each neuron in each cluster - staggering
%locations
for i = 1:clusters
    for j = 1:n
        if i <= clusters/2
            x = i;
            y_base = (clusters/2) - abs((clusters/4) - i);
        else
            x = clusters - i;
            y_base = abs((3*clusters/4) - i);
        end
        cluster_plot_indices(i,j,1) = x + 0.1*randn();
        cluster_plot_indices(i,j,2) = y_base + 0.1*randn();
    end
end

%Next find overlapping neuron scatter positions to draw connecting lines
cluster_plot_connections = [];
for i = 1:n
    clusters_in = find(cluster_mat(:,i));
    if length(clusters_in) > 2
        possible_pairs = nchoosek(clusters_in,2);
        for j = 1:length(possible_pairs)
            index_pair = zeros(2,2);
            index_pair(1,:) = cluster_plot_indices(possible_pairs(j,1),i,:);
            index_pair(2,:) = cluster_plot_indices(possible_pairs(j,2),i,:);
            cluster_plot_connections = cat(3,cluster_plot_connections,index_pair); %#ok<AGROW>
        end    
    elseif length(clusters_in) == 2
        index_pair = zeros(2,2);
        index_pair(1,:) = cluster_plot_indices(clusters_in(1),i,:);
        index_pair(2,:) = cluster_plot_indices(clusters_in(2),i,:);
    end
end
[~, ~, cluster_connections_l] = size(cluster_plot_connections);

%Set colors for each cluster
color_options = jet(clusters);

%Plot
f = figure;
leg_items = [];
leg_names = [];
hold on
for i = 1:clusters
    scat = scatter(cluster_plot_indices(i,:,1),cluster_plot_indices(i,:,2),[],color_options(i,:),'filled');
    leg_items(end+1) = scat;
    leg_names = [leg_names, strcat('Cluster #',string(i))]; %#ok<AGROW>
end
for i = 1:cluster_connections_l
    line(cluster_plot_connections(:,1,i),cluster_plot_connections(:,2,i))
end
legend(leg_items,leg_names)
title('Network Diagram','FontSize',16)
set(gca,'xtick',[],'ytick',[])


%Plot individual cluster connections
f = figure;
hold on
for i = 1:clusters
    within_cluster_connections = [];
    %First find and store intra-cluster connections
    neur_ind = find(cluster_mat(i,:));
    for j = neur_ind %row
        for k = neur_ind %column
            if conns(j,k) == 1
                index_pair = zeros(2,2);
                index_pair(1,:) = cluster_plot_indices(i,j,:);
                index_pair(2,:) = cluster_plot_indices(i,k,:);
                within_cluster_connections = cat(3,within_cluster_connections,index_pair); %#ok<AGROW>
            end    
        end
    end
    [~, ~, cluster_connections_l] = size(within_cluster_connections);
    %Plot
    scat = scatter(cluster_plot_indices(i,:,1),cluster_plot_indices(i,:,2),[],color_options(i,:),'filled');
    for i = 1:cluster_connections_l
        line(within_cluster_connections(:,1,i),within_cluster_connections(:,2,i))
    end
end 
title('Intra-Cluster Connectivity','FontSize',16)
set(gca,'xtick',[],'ytick',[])

