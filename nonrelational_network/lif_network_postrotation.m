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
tau_syn_E = 10*10^(-3); %AMPA/NMDA synaptic decay time constant (s)
tau_syn_I = 5*10^(-3); %GABA synaptic decay time constant (s)
tau_sra = 30*10^(-3); %spike rate adaptation time constant (s)
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
del_G_syn_E = 8*10^(-9); %synaptic conductance step following spike (S)
del_G_syn_I = 8*10^(-9); %synaptic conductance step following spike (S)
%___________________________
del_G_sra = 200*10^(-9); %spike rate adaptation conductance step following spike %ranges from 1-200 *10^(-9) (S)
connectivity_gain = 1; %amount to increase or decrease connectivity by with each spike (more at the range of 1.002-1.005)
IEI = 0.1; %inter-event-interval (s) the elapsed time between spikes to count separate events
IES = ceil(IEI/dt); %inter-event-steps = the number of steps to elapse between spikes

%How spikes are initiated:
%'cluster' sets a cluster to threshold;
%'current' means spikes depend on an input current of I_in; 
%'neuron' sets a random selection of 2% of neurons to threshold
type = 'current'; %'neuron'; %'cluster';

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
I_coeff = 2.7; %set to 0 for no input current
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

%Save parameters to a structure and to computer
w = whos;
parameters = struct;
for a = 1:length(w) 
    parameters.(w(a).name) = eval(w(a).name); 
end
clear w a
save(strcat(save_path,'/parameters.mat'),'parameters')

%% Create networks and test spike progressions
%Runs through a series of different random number generator seeds to change
%the network connectivity and setup, and then automatically outputs
%sequence data to a folder

%If uploading a parameter file, uncomment the next line
load(strcat(save_path,'/parameters.mat'))

for i = 1:10 %how many different network structures to test
    rng(i) %set random number generator for network structure
    
    %CREATE NETWORK SAVE PATH
    net_save_path = strcat(save_path,'/network_',string(i));
    if ~isfolder(net_save_path)
        mkdir(net_save_path);
    end
    
    %SET UP NETWORK
    [cluster_mat, conns] = create_clusters(parameters.n, ...
        parameters.clusters, parameters.cluster_n, parameters.cluster_prob);
    conns_copy = conns; %just a copy of the connections to maintain for reset runs if there's "plasticity"
    %Randomize excitatory and inhibitory connection strengths based on selected
    %probability.
    all_indices = [1:parameters.n];
    I_indices = randi(parameters.n,[parameters.n_I,1]); %indices of inhibitory neurons
    all_indices(I_indices) = 0;
    E_indices = find(all_indices)'; %indices of excitatory neurons
    n_EE = sum(conns(E_indices,E_indices),'all'); %number of E-E connections
    n_EI = sum(conns(E_indices,I_indices),'all'); %number of E-I connections
    n_II = sum(conns(I_indices,I_indices),'all'); %number of I-I connections
    n_IE = sum(conns(I_indices,E_indices),'all'); %number of I-E connections
    clear all_indices
    
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
    network_spike_sequences = struct;
    network_cluster_sequences = struct; 
    
    for j = 1:parameters.test_val_max
        seed = j;
        
        %Create Storage Variables
        I_syn = zeros(parameters.n,parameters.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
        %synaptic conductances for each neuron at each timestep
        G_syn_I = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic inhibitory (S)
        G_syn_E = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic excitatory (S)
        V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
        G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep
        
        %Run model
        [V_m, G_sra, G_syn_I, G_syn_E, I_syn] = lif_sra_calculator_postrotation(...
            parameters, seed, cluster_mat, conns, I_indices, E_indices, ...
            I_syn, G_syn_I, G_syn_E, V_m, G_sra);
        network_var(j).V_m = V_m;
        network_var(j).G_sra = G_sra;
        network_var(j).G_syn_I = G_syn_I;
        network_var(j).G_syn_E = G_syn_E;
        network_var(j).I_syn = I_syn;
        
        %Find spike profile
        spikes_V_m = V_m >= parameters.V_th;
        [spikes_x,spikes_t] = find(spikes_V_m);
        max_time = max(spikes_t);
        
        %Find maximum firing rate + average maximum firing rates of neurons
        all_fr = sum(spikes_V_m,2)/parameters.t_max;
        max_fr = max(all_fr);
        avg_fr = mean(all_fr);
        
        %USE FR AS AN INDICATOR OF GOOD SEQUENCES
        %In hippocampus, Jadhav lab uses 3 Hz as a place-cell cutoff
        if and(max_fr <= 3, max_fr >= 1) && and(avg_fr <= 2, avg_fr >= 1/(parameters.t_max+1))
            %Find event times
            events = []; 
            last_start = spikes_t(1);
            last_time = spikes_t(1);
            spike_count = 1;
            for t_i = 2:length(spikes_t)
                s_i = spikes_t(t_i);
                if s_i - last_time <= parameters.IES
                    last_time = s_i;
                    spike_count = spike_count + 1;
                else
                    if (last_start ~= last_time) && (spike_count > 0.05*parameters.n) %weed out events of single spikes
                        events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last range of spikes to the events vector
                    end
                    last_start = s_i;
                    last_time = s_i;
                    spike_count = 1;
                end
            end
            if last_start ~= last_time %if the last event is a single spike, we don't care about it
                events = [events; [last_start, last_time]]; %#ok<AGROW> %add the last interval
            end
            [num_events,~] = size(events);
            display(events)
            %save to both structures
            network_spike_sequences(j).events = events;
            network_cluster_sequences(j).events = events;

            %Find spike sequences
            for e_i = 1:num_events
                %store spike orders for each event
                event_spikes = spikes_V_m(:,events(e_i,1):events(e_i,2));
                [e_spikes_x, ~] = find(event_spikes);
                spike_order = unique(e_spikes_x,'stable');
                network_spike_sequences(j).spike_order.(strcat('sequence_',string(e_i))) = spike_order;
                %store ranks for each neuron
                ranks_vec = zeros(1,parameters.n);
                for k = 1:length(spike_order)
                    n_ind = spike_order(k);
                    ranks_vec(1,n_ind) = k;
                end
                network_spike_sequences(j).spike_ranks.(strcat('sequence_',string(e_i))) = ranks_vec;
                %store nonspiking neurons
                nonspiking_neurons = isnan(ranks_vec./ranks_vec);
                network_spike_sequences(j).nonspiking_neurons.(strcat('sequence_',string(e_i))) = nonspiking_neurons;
            end
            clear e_i event_spikes e_spikes_x spike_order ranks_vec k n_ind nonspiking_neurons
            
            %Visualize re-ordered spike sequences
            f = figure;
            axes = [];
            for e_i = 1:num_events
                spike_order = network_spike_sequences(j).spike_order.(strcat('sequence_',string(e_i)));
                sub_spikes_V_m = spikes_V_m(:,events(e_i,1):events(e_i,2));
                reordered_spikes = zeros(size(sub_spikes_V_m));
                [~,event_length] = size(reordered_spikes);
                for s_i = 1:length(spike_order)
                    reordered_spikes(s_i,:) = sub_spikes_V_m(spike_order(s_i),:);
                end 
                ax = subplot(1,num_events,e_i);
                axes = [axes, ax];
                imagesc(reordered_spikes)
                xticks(round(linspace(1,event_length,20))) %20 ticks will be displayed
                xt = get(gca,'XTick');
                xtlbl = round(linspace(events(e_i,1)*parameters.dt,events(e_i,2)*parameters.dt,numel(xt)),2);
                colormap(flip(gray))
                xlabel('Time (s)','FontSize',16)
                ylabel('Reordered Neuron Number','FontSize',16)
                title(strcat('Event #',string(e_i)))
                set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
            end
            if strcmp(type,'cluster')
                sgtitle(strcat('Spiking Behavior: cluster = ',string(j)),'FontSize',16)
            else
                sgtitle('Spiking Behavior','FontSize',16)
            end
            %linkaxes(axes)
            savefig(f,strcat(net_save_path,'/',type,'_',string(j),'firing_sequence.fig'))
            saveas(f,strcat(net_save_path,'/',type,'_',string(j),'firing_sequence.jpg'))
            close(f)
            clear e_i spike_order reordered_spikes event_length s_i ax ...
                axes xt xtlbl
            
            %Find cluster sequence per event by moving bin
            bin_width = 5*10^(-3); %5 ms bin
            bin_size = ceil(bin_width/parameters.dt); %number of timesteps to use in a bin
            for e_i = 1:num_events
                cluster_spikes = cluster_mat*spikes_V_m(:,events(e_i,1):events(e_i,2));
                cluster_mov_sum = movsum(cluster_spikes',bin_size)';
                normalized_cluster_spikes = cluster_spikes ./ sum(cluster_spikes,1);
                normalized_cluster_spikes(isnan(normalized_cluster_spikes)) = 0;
                normalized_cluster_mov_sum = cluster_mov_sum ./ sum(cluster_mov_sum,1);
                normalized_cluster_mov_sum(isnan(normalized_cluster_mov_sum)) = 0;
                network_cluster_sequences(j).clusters.(strcat('sequence_',string(e_i))) = cluster_spikes;
                network_cluster_sequences(j).movsum.(strcat('sequence_',string(e_i))) = cluster_mov_sum;
                network_cluster_sequences(j).normalized_clusters.(strcat('sequence_',string(e_i))) = normalized_cluster_spikes;
                network_cluster_sequences(j).normalized_cluster_mov_sum.(strcat('sequence_',string(e_i))) = normalized_cluster_mov_sum;
            end
            clear bin_width bin_size cluster_spikes cluster_mov_sum e_i
        end
    end
    
    %SAVE NETWORK DATA
    save(strcat(net_save_path,'/network_var.mat'),'network_var','-v7.3')
    save(strcat(net_save_path,'/network_spike_sequences.mat'),'network_spike_sequences','-v7.3')
    save(strcat(net_save_path,'/network_cluster_sequences.mat'),'network_cluster_sequences','-v7.3')
end

%% Calculate trajectory similarity
%This code block does not require any of the earlier code blocks to run. It
%simply requires stored data from prior runs of earlier code blocks.

%Select and load data to analyze
net_save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Network Save Folder'); %Have user input where they'd like the output stored
slashes = find(net_save_path == '/');
save_path = net_save_path(1:slashes(end));
load(strcat(save_path,'/parameters.mat'))
load(strcat(net_save_path,'/network_spike_sequences.mat'))

[~,inits] = size(network_spike_sequences);

[neuron_ranks, rho_unique, ranks, rho_nonunique, ranks_mod] = calculate_trajectory_similarity(n, ...
    inits, network_spike_sequences);

ranks(isnan(ranks)) = 0;
ranks_mod(isnan(ranks_mod)) = 0;

average_rank = mean(ranks(ranks~=1),'all');
std_rank = std(ranks(ranks~=1),[],'all');
average_rank_mod = mean(ranks_mod(ranks_mod~=1),'all');
std_rank_mod = std(ranks_mod(ranks_mod~=1),[],'all');

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

