%% Initialize Variables

%close all
clear
rng(1) %same randomizations each time

%Neuron and cluster counts
n = 200; %number of neurons
clusters = round(n/20); %number of clusters of neurons (for small n round(n/5), for large n round(n/20)) 
cluster_n = round(n/5); %number of neurons in a cluster (for small n round(n/3), for large n round(n/5)) 

%Interaction constants
t_max = 2; %maximum amount of time (s)
dt = 0.1*10^(-3); %timestep (s)
t_steps = t_max/dt;
tau_syn_E = 20*10^(-3); %AMPA/NMDA synaptic decay time constant (s)
tau_syn_I = 5*10^(-3); %GABA synaptic decay time constant (s)
tau_sra = 30*10^(-3); %spike rate adaptation time constant (s)
E_K = -80*10^(-3); %potassium reversal potential (V)
E_L = -70*10^(-3); %leak reversal potential (V)
G_L = 25*10^(-9); %leak conductance (S)
C_m = 0.5*10^(-9); %total membrane capacitance (F)
V_th = -50*10^(-3); %threshold membrane potential (V)
V_reset = -70*10^(-3); %reset membrane potential
V_syn_E = 0; %synaptic reversal potential (excitatory)
V_syn_I = -70*10^(-3); %synaptic reversal potential (inhibitory)
E_syn_E = V_syn_E*ones(n,1); %vector of the synaptic reversal potential for excitatory connections
E_syn_I = V_syn_I*ones(n,1); %vector of the synaptic reversal potential for inhibitory connections
%______Split del_G_syn______
del_G_syn_E = 15*10^(-9); %excitatory synaptic conductance step following spike
del_G_syn_I = 10*10^(-9); %inhibitory synaptic conductance step following spike
%___________________________
del_G_sra = 200*10^(-9); %spike rate adaptation conductance step following spike
connectivity_gain = 1.005; %amount to increase or decrease connectivity by with each spike (more at the range of 1.002-1.005)
%Added theta wave input
x_in = [0:dt:t_max];
I_coeff = 5.1*10^(-10);
I_Hz = 8;
I_theta = I_coeff*(0.5 + 0.5*sin(I_Hz*2*pi*x_in)); %zeros(1,t_steps+1); %Approximately theta wave input current

%Calculate connection probabilites
conn_prob = 0.1; %set a total desired connection probability
npairs = n*(n-1); %total number of possible neuron connections
nclusterpairs = cluster_n*(cluster_n - 1)*clusters; %total number of possible intra-cluster connections
cluster_prob = min(conn_prob*nclusterpairs/npairs,1); %intra-cluster connection probability

%Create clusters (either totally random or spacially clustered)
cluster_mat = zeros(clusters,n);
for i = 1:clusters %set clusters
    ord = randperm(n,cluster_n); %random connectivity
    cluster_mat(i,ord) = 1; %note which neurons are in a cluster
end
[cluster_number,neuron_number] = find(cluster_mat);

clear ord i
conn_bin = (cluster_mat'*cluster_mat) > 0; %binary matrix of connections
for i = 1:n
    conn_bin(i,i) = 0; %remove self connectivity
end
syn_count = sum(conn_bin,2); %number of connections each neuron has
non_connected = find(syn_count == 0); %check if any are not in a cluster
% for i = non_connected %add the non-connected back in
%     rand_clust = randi(cluster_n);
%     cluster_mat(rand_clust,i) = 1;
% end
%conns = cluster_mat'*cluster_mat; %let all connections exist
conns = zeros([n,n]);
for i = 1:clusters %now set the connectivity matrix using probabilities
    ord = find(cluster_mat(i,:));
    conns(ord,ord) = conns(ord,ord) + 0.2*(rand(length(ord)) < cluster_prob); %set connection probabilities (all positive)
end
for i = 1:n
    conns(i,i) = 0; %remove self connectivity
end
clear i rand_clust ord 
conn_bin = conns > 0; %binary matrix of locations of connectivity
conns_copy = conns; %just a copy of the connections to maintain for reset runs if there's "plasticity"
syn_count_in = sum(conn_bin,2); %number of connections each neuron receives
syn_count_out = sum(conn_bin,1)'; %number of connections each neuron makes
overall_connected = (syn_count_in > 0) + (syn_count_out > 0);

%Randomize excitatory and inhibitory connection strengths based on selected
%probability.
p_E = 0.8; %probability of an excitatory neuron
p_I = 1 - p_E; %probability of an inhibitory neuron
n_I = round(p_I*n); %number of inhibitory neurons
all_indices = [1:n];
I_indices = randi(n,[n_I,1]); %indices of inhibitory neurons
all_indices(I_indices) = 0;
E_indices = find(all_indices)'; %indices of excitatory neurons
clear all_indices n_I
n_EE = sum(conn_bin(E_indices,E_indices),'all'); %number of E-E connections
n_EI = sum(conn_bin(E_indices,I_indices),'all'); %number of E-I connections
n_II = sum(conn_bin(I_indices,I_indices),'all'); %number of I-I connections
n_IE = sum(conn_bin(I_indices,E_indices),'all'); %number of I-E connections

%Create storage variables
I_syn = zeros(n,t_steps+1); %synaptic current emitted by each neuron at each timestep
%synaptic conductances for each neuron at each timestep
G_syn_I = zeros(n,t_steps+1); %conductance for presynaptic inhibitory
G_syn_E = zeros(n,t_steps+1); %conductance for presynaptic excitatory
V_m = zeros(n,t_steps+1); %membrane potential for each neuron at each timestep
G_sra = zeros(n,t_steps+1); %refractory conductance for each neuron at each timestep

V_m(:,1) = V_reset + randn([n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise

clear i neur_start

%% Run Model Using Functions

%Set the type of firing.
%   'neuron' randomly selects 'cluster' number of neurons to start firing
%   'cluster' randomly selects a cluster and sets all of the cluster
%       neurons to start firing
type = 'cluster'; %'neuron'; 

%Set whether neurons reset to baseline at a changepoint
%   1 means reset to baseline of no firing for all neurons except those
%   randomly selected. 0 means do not reset and just use the firing rate
%   from the previous step.
reset = 1;

%Set the number of times to introduct a changepoint
change_num = 0;

%Set seed for the random number generator (for first spiking neurons)
seed = 5;

%Run with changepoints representing complete activity stops and start of
%new environment stimulation
if change_num > 0
    change_pts = round(linspace(1,t_steps,change_num+1));
    for i = 2:change_num+1
        t_start = change_pts(i-1);
        if reset == 1
            t_end = change_pts(i) - 1; %this ensures the next timestep is all 0s for other variables
            conns = conns_copy; %reset connections
        else
            t_end = change_pts(i);
        end
        [V_m, G_sra, G_syn_I, G_syn_E, I_syn, conns] = lif_sra_calculator(n, ...
            clusters,cluster_mat,conns,conn_bin,V_m,V_reset,V_th,G_sra,del_G_sra,G_syn_I, ...
            G_syn_E, del_G_syn_E,del_G_syn_I,I_syn,E_syn_E,E_syn_I,E_K,E_L,G_L,C_m,t_start,t_end,dt, ...
            tau_sra,tau_syn_E,tau_syn_I,connectivity_gain,type,reset,I_indices,E_indices,I_theta,seed);
    end
else
    t_start = 1;
    t_end = t_steps;
    change_pts = [t_start,t_end];
    [V_m, G_sra, G_syn_I, G_syn_E, I_syn, conns] = lif_sra_calculator(n, ...
        clusters,cluster_mat,conns,conn_bin,V_m,V_reset,V_th,G_sra,del_G_sra,G_syn_I, ...
        G_syn_E, del_G_syn_E,del_G_syn_I,I_syn,E_syn_E,E_syn_I,E_K,E_L,G_L,C_m,t_start,t_end,dt, ...
        tau_sra,tau_syn_E,tau_syn_I,connectivity_gain,type,reset,I_indices,E_indices,I_theta,seed);
end

%Plot membrane potentials of all neurons over time
%As a line plot
% figure;
% plot(V_m')
% hold on
% for i = 1:length(change_pts)
%     xline(change_pts(i),':','LineWidth',2)
% end
% clear i
% xlabel('Time Step','FontSize',16)
% ylabel('Membrane Potential (V)','FontSize',16)
% title('Membrane Potential of All Neurons')

%This section should be uncommented only after the "Re-order" code block is run
%and tests of different seeds are being run in this section using the
%original seed's spike ordering.
% reordered_spikes_2 = zeros(size(spikes_V_m));
% for i = 1:n
%     reordered_spikes_2(i,:) = spikes_V_m(spike_order(i),:);
% end
% spike_time_indices_2 = sum(reordered_spikes_2,1) > 0;
% spikes_collapsed_2 = reordered_spikes_2(:,spike_time_indices);
% figure;
% imagesc(spikes_collapsed_2)
% colormap(flip(gray))
% xlabel('Time Step','FontSize',16)
% ylabel('Reordered Neuron Number','FontSize',16)
% %seed 1, cluster 5; seed 3, cluster 6; seed 5, cluster 3
% title('Spiking Behavior (Collapsed): seed = 5, cluster = 3','FontSize',16)

%% Re-order neurons based on spike time and show plot
%If testing how reordering for one seed looks for other seeds, after running
%this section, uncomment last bit of code in the above code block and re-run
%it for different seeds.

%Spike Times
spikes_V_m = V_m >= V_th;
[spikes_x,spikes_t] = find(spikes_V_m);

%Store the order of spikes
spike_order = [];
for i = 1:length(spikes_x) %run through spiking neuron order
    if isempty(find(spike_order == spikes_x(i),1))
        spike_order(end+1) = spikes_x(i); %#ok<*SAGROW>
    end
end
for i = 1:n %make sure all indices are accounted for (non-spiking neurons)
    if isempty(find(spike_order == i,1))
        spike_order(end+1) = i;
    end
end

%Show in imagesc the re-ordered neurons' spiking
reordered_spikes = zeros(size(spikes_V_m));
reordered_V_m = zeros(size(V_m));
for i = 1:n
    reordered_spikes(i,:) = spikes_V_m(spike_order(i),:);
    reordered_V_m(i,:) = V_m(spike_order(i),:);
end
figure;
imagesc(reordered_spikes)
% colormap(flip(gray))
% title('Reordered Spiking Behavior','FontSize',16)
title(strcat('Spiking Behavior: seed = ',string(seed),' I_{theta} = ', string(I_coeff)),'FontSize',16)
xlabel('Time Step','FontSize',16)
ylabel('Reordered Neuron Number','FontSize',16)
% xlim([450,1000])

%to visualize breaks between theta waves in collapsed view
theta_breaks = find(round(I_theta,13) == round(I_coeff/2,13));
reord_spikes_lined = reordered_spikes;
reord_spikes_lined(:,theta_breaks) = reordered_spikes(:,theta_breaks) + 0.25;

% %Collapsed reordered spikes
spike_time_indices = sum(reord_spikes_lined,1) > 0; %sum(reordered_spikes,1) > 0;
spikes_collapsed = reord_spikes_lined(:,spike_time_indices); %reordered_spikes(:,spike_time_indices);
figure;
imagesc(spikes_collapsed)
colormap(flip(gray))
%seed 1, cluster 5; seed 3, cluster 6; seed 5, cluster 3
title(strcat('Reordered Spiking Behavior (Collapsed): seed = ',string(seed),' I_{theta} = ', string(I_coeff)),'FontSize',16)
% title('Reordered Spiking Behavior (Collapsed)')

%Visual of membrane potential of reordered neurons over time
figure;
imagesc(reordered_V_m)

%% Run Model and store ordering for all cluster initiations
%Prior to running this section, must update lif_sra_calculator.m code to
%have i = seed in line 30.
%Make sure to save all_cluster_data and neuron_reorder after running for
%re-use in the future.

all_cluster_data = struct; %store all statistics for each cluster as an initiator

type = 'cluster';
reset = 1;
change_num = 0;
neuron_reorder = zeros(n,clusters);

for seed = 1:clusters
    %Run model
    t_start = 1;
    t_end = t_steps;
    [V_m, G_sra, G_syn_I, G_syn_E, I_syn, conns] = lif_sra_calculator(n, ...
        clusters,cluster_mat,conns,conn_bin,V_m,V_reset,V_th,G_sra,del_G_sra,G_syn_I, ...
        G_syn_E, del_G_syn_E,del_G_syn_I,I_syn,E_syn_E,E_syn_I,E_K,E_L,G_L,C_m,t_start,t_end,dt, ...
        tau_sra,tau_syn_E,tau_syn_I,connectivity_gain,type,reset,I_indices,E_indices,I_theta,seed);
    
    %Collect basic data
    spikes_V_m = V_m >= V_th;
    [spikes_x,spikes_t] = find(spikes_V_m);
    spike_count = sum(spikes_V_m,2);
    
    %Find new spike ordering
    spike_order = [];
    for i = 1:length(spikes_x) %run through spiking neuron order
        if isempty(find(spike_order == spikes_x(i),1))
            spike_order(end+1) = spikes_x(i); %#ok<*SAGROW>
        end
    end
    non_spiking = [];
    for i = 1:n %make sure all indices are accounted for (non-spiking neurons)
        if isempty(find(spike_order == i,1))
            spike_order(end+1) = i;
            non_spiking(end+1) = i;
        end
    end
    
    %store data
    all_cluster_data(seed).V_m = V_m;
    all_cluster_data(seed).G_sra = G_sra; 
    all_cluster_data(seed).G_syn_I = G_syn_I;
    all_cluster_data(seed).G_syn_E = G_syn_E;
    all_cluster_data(seed).I_syn = I_syn;
    all_cluster_data(seed).conns = conns;
    all_cluster_data(seed).spike_count = spike_count;
    all_cluster_data(seed).spikes = [spikes_x,spikes_t];
    all_cluster_data(seed).neuron_reorder = spike_order;
    all_cluster_data(seed).non_spiking = non_spiking;
    neuron_reorder(:,seed) = spike_order';
end

clear non_spiking i spike_order t_start t_end seed
%% Analyze correlation between different seed spiking orders
%This section needs variables from above section (all_cluster_data and 
%neuron_reorder) and the variable initiation section to run.

cluster_pairs = nchoosek([1:clusters],2);

%calculate ranks
neuron_ranks = zeros(n,clusters);
for i = 1:clusters
    for j = 1:n
        ind = find(neuron_reorder(:,i) == j);
        neuron_ranks(j,i) = ind;
    end    
end

%Using Spearman's rank correlation on the unique integer values we get
rho_unique = corr(neuron_ranks,'Type','Spearman');

%If all n ranks are distinct integers, can use the formula: (wikipedia)
%r_s = 1 - (6*sum(d_i^2))/(n*(n^2 - 1))
%d_i = rg(X_i) - rg(Y_i)
ranks = zeros(clusters,clusters);
for i = 1:clusters
    X_i = neuron_ranks(:,i);
    for j = 1:clusters
        Y_i = neuron_ranks(:,j);
        d = zeros(1,n);
        for k = 1:n
            ix_i = find(X_i == k);
            iy_i = find(Y_i == k);
            d(1,k) = ix_i - iy_i;
        end
        ranks(i,j) = 1 - (6*sum(d.^2))/(n*(n^2 - 1));
        ranks(j,i) = ranks(cluster_pairs(i,1),cluster_pairs(i,2));
        if i == j
            ranks(i,j) = 1;
        end
    end
end
clear i X_i Y_i d j ix_i iy_i k
average_rank = mean(ranks(ranks~=1),'all');
std_rank = std(ranks(ranks~=1),[],'all');

%Above rank calculation modified to exclude non-spiking neurons
ranks_mod = zeros(clusters,clusters);
for i = 1:clusters
    X_i = neuron_ranks(:,i);
    X_nonspike = all_cluster_data(i).non_spiking;
    for j = 1:clusters
        Y_i = neuron_ranks(:,j);
        Y_nonspike = all_cluster_data(j).non_spiking;
        nonspiking = unique([X_nonspike,Y_nonspike]);
        %To ensure we are only comparing neurons that fire in both simulations,
        %we have to remove nonspiking neurons from both sets.
        X_i = setdiff(X_i, nonspiking','stable');
        Y_i = setdiff(Y_i, nonspiking','stable');
        new_n = (length(X_i)+length(Y_i))/2;
        d = zeros(1,n);
        for k = 1:n
            ix_i = find(X_i == k);
            iy_i = find(Y_i == k);
            if ~isempty(ix_i) & ~isempty(iy_i)
                d(1,k) = ix_i - iy_i;
            end
        end
        ranks_mod(i,j) = 1 - (6*sum(d.^2))/(new_n*(new_n^2 - 1));
        ranks_mod(j,i) = ranks_mod(i,j);
        if i == j
            ranks_mod(i,j) = 1;
        end
    end
end
clear i X_i Y_i d j new_n k ix_i iy_i X_nonspike Y_nonspike
average_rank_mod = mean(ranks_mod(ranks_mod~=1),'all');
std_rank_mod = std(ranks_mod(ranks_mod~=1),[],'all');


%From Vas et al. Matching Index calculation
matching_index = zeros(clusters,clusters);
for clust_pair = 1:length(cluster_pairs) %for each pair of clusters calculate MI
    match = 0;
    non_match = 0;
    c_1 = cluster_pairs(clust_pair,1);
    c_2 = cluster_pairs(clust_pair,2);
    for i = 1:n %first neuron index
        [c1_i,~] = find(neuron_ranks(:,c_1) == i);
        [c2_i,~] = find(neuron_ranks(:,c_2) == i);
        for j = i:n %second neuron index
            [c1_j,~] = find(neuron_ranks(:,c_1) == j);
            [c2_j,~] = find(neuron_ranks(:,c_2) == j);
            %Calculate the logical values of whether i is after j for both
            %cluster initiation regimes
            g1_i = c1_i > c1_j;
            g2_i = c2_i > c2_j;
            if g1_i == g2_i %if both have the same ordering, they will be equal (either both 1 or both 0)
                match = match + 1;
            else
                non_match = non_match + 1;
            end
        end
    end 
    matching_index(c_1,c_2) = (match - non_match)/ (match + non_match);
    matching_index(c_2,c_1) = matching_index(c_1,c_2);
end
clear clust_pair match non_match c_1 c_2 i j c1_i c2_i c1_j c2_j
average_matching_index = mean(matching_index(matching_index~=1),'all');
std_matching_index = std(matching_index(matching_index~=1),[],'all');

%From Vas et al. Matching Index calculation (MODIFIED - does not count
%neurons that both do not spike)
matching_index_mod = zeros(clusters,clusters);
for clust_pair = 1:length(cluster_pairs) %for each pair of clusters calculate MI
    match = 0;
    non_match = 0;
    c_1 = cluster_pairs(clust_pair,1);
    spikes_x_1 = all_cluster_data(c_1).spikes(:,1);
    c_2 = cluster_pairs(clust_pair,2);
    spikes_x_2 = all_cluster_data(c_2).spikes(:,1);
    for i = 1:n %first neuron index
        [c1_i,~] = find(spikes_x_1 == i,1);
        [c2_i,~] = find(spikes_x_2 == i,1);
        for j = i:n %second neuron index
            [c1_j,~] = find(spikes_x_1 == j,1);
            [c2_j,~] = find(spikes_x_2 == j,1);
            %Calculate the logical values of whether i is after j for both
            %cluster initiation regimes
            if isempty(c1_i) && ~isempty(c1_j)
                g1_i = 0;
                if isempty(c2_i) && ~isempty(c2_j)
                    g2_i = 0;
                    match = match+1;
                else
                    non_match = non_match + 1;
                end
            elseif isempty(c1_j) && ~isempty(c1_i)
                g1_i = 1;
                if isempty(c2_j) && ~isempty(c2_i)
                    g2_i = 1;
                    match = match+1;
                else
                    non_match = non_match + 1;
                end
            elseif ~isempty(c1_i) && ~isempty(c1_j)
                g1_i = c1_i > c1_j;
                if ~isempty(c2_i) && ~isempty(c2_j)
                    g2_i = c2_i > c2_j;
                    if g1_i == g2_i %if both have the same ordering, they will be equal (either both 1 or both 0)
                        match = match + 1;
                    else
                        non_match = non_match + 1;
                    end
                else
                    non_match = non_match + 1;
                end
            end
        end
    end 
    matching_index_mod(c_1,c_2) = (match - non_match)/ (match + non_match);
    matching_index_mod(c_2,c_1) = matching_index_mod(c_1,c_2);
end
clear clust_pair match non_match c_1 c_2 i j c1_i c2_i c1_j c2_j
average_matching_index_mod = mean(matching_index_mod(matching_index_mod~=1),'all');
std_matching_index_mod = std(matching_index_mod(matching_index_mod~=1),[],'all');

%% Calculate Matching Index Distributions for Shuffles

%From Vas et al. Matching Index calculation with shuffled data
matching_index = zeros(clusters,clusters);
for clust_pair = 1:length(cluster_pairs) %for each pair of clusters calculate MI
    match = 0;
    non_match = 0;
    c_1 = cluster_pairs(clust_pair,1);
    c_2 = cluster_pairs(clust_pair,2);
    for i = 1:n %first neuron index
        [c1_i,~] = find(neuron_ranks(:,c_1) == i);
        [c2_i,~] = find(neuron_ranks(:,c_2) == i);
        for j = i:n %second neuron index
            [c1_j,~] = find(neuron_ranks(:,c_1) == j);
            [c2_j,~] = find(neuron_ranks(:,c_2) == j);
            %Calculate the logical values of whether i is after j for both
            %cluster initiation regimes
            g1_i = c1_i > c1_j;
            g2_i = c2_i > c2_j;
            if g1_i == g2_i %if both have the same ordering, they will be equal (either both 1 or both 0)
                match = match + 1;
            else
                non_match = non_match + 1;
            end
        end
    end 
    matching_index(c_1,c_2) = (match - non_match)/ (match + non_match);
    matching_index(c_2,c_1) = matching_index(c_1,c_2);
end
clear clust_pair match non_match c_1 c_2 i j c1_i c2_i c1_j c2_j
average_matching_index = mean(matching_index(matching_index~=1),'all');
std_matching_index = std(matching_index(matching_index~=1),[],'all');

%From Vas et al. Matching Index calculation (MODIFIED - does not count
%neurons that both do not spike)
matching_index_mod = zeros(clusters,clusters);
for clust_pair = 1:length(cluster_pairs) %for each pair of clusters calculate MI
    match = 0;
    non_match = 0;
    c_1 = cluster_pairs(clust_pair,1);
    spikes_x_1 = all_cluster_data(c_1).spikes(:,1);
    c_2 = cluster_pairs(clust_pair,2);
    spikes_x_2 = all_cluster_data(c_2).spikes(:,1); 
    spikes_x_2 = spikes_x_2(randperm(length(spikes_x_2))); %shuffled
    for i = 1:n %first neuron index
        [c1_i,~] = find(spikes_x_1 == i,1);
        [c2_i,~] = find(spikes_x_2 == i,1);
        for j = i:n %second neuron index
            [c1_j,~] = find(spikes_x_1 == j,1);
            [c2_j,~] = find(spikes_x_2 == j,1);
            %Calculate the logical values of whether i is after j for both
            %cluster initiation regimes
            if isempty(c1_i) && ~isempty(c1_j)
                g1_i = 0;
                if isempty(c2_i) && ~isempty(c2_j)
                    g2_i = 0;
                    match = match+1;
                else
                    non_match = non_match + 1;
                end
            elseif isempty(c1_j) && ~isempty(c1_i)
                g1_i = 1;
                if isempty(c2_j) && ~isempty(c2_i)
                    g2_i = 1;
                    match = match+1;
                else
                    non_match = non_match + 1;
                end
            elseif ~isempty(c1_i) && ~isempty(c1_j)
                g1_i = c1_i > c1_j;
                if ~isempty(c2_i) && ~isempty(c2_j)
                    g2_i = c2_i > c2_j;
                    if g1_i == g2_i %if both have the same ordering, they will be equal (either both 1 or both 0)
                        match = match + 1;
                    else
                        non_match = non_match + 1;
                    end
                else
                    non_match = non_match + 1;
                end
            end
        end
    end 
    matching_index_mod(c_1,c_2) = (match - non_match)/ (match + non_match);
    matching_index_mod(c_2,c_1) = matching_index_mod(c_1,c_2);
end
clear clust_pair match non_match c_1 c_2 i j c1_i c2_i c1_j c2_j
average_matching_index_mod = mean(matching_index_mod(matching_index_mod~=1),'all');
std_matching_index_mod = std(matching_index_mod(matching_index_mod~=1),[],'all');

%% Analyze similarity between different seeds (general)
%needs variables from above section (all_cluster_data and neuron_reorder)

cluster_pairs = nchoosek([1:clusters],2);

%spiking similarity: calculate how similar the spiking neurons are from one
%simulation to another (fraction)
spiking_similarity = zeros(clusters,clusters);
for i = 1:clusters
    spike_i = unique(all_cluster_data(i).spikes(:,1));
    for j = 1:clusters
        spike_j = unique(all_cluster_data(j).spikes(:,1));
        matched_spiking = setdiff(spike_i, spike_j);
        spiking_similarity(i,j) = 1 - length(matched_spiking)/(0.5*(length(spike_i) + length(spike_j)));
    end
end
clear i spike_i j spike_j matched_spiking
avg_spike_sim = mean(spiking_similarity(spiking_similarity ~= 1),'all');
std_spike_sim = std(spiking_similarity(spiking_similarity ~= 1),[],'all');

%non-spiking similarity: calculate how similar the non-spiking neurons are
%from one simulation to another (fraction)
non_spiking_similarity = zeros(clusters,clusters);
for i = 1:clusters
    non_spike_i = all_cluster_data(i).non_spiking;
    for j = 1:clusters
        non_spike_j = all_cluster_data(j).non_spiking;
        matched_nonspiking = setdiff(non_spike_i, non_spike_j);
        non_spiking_similarity(i,j) = 1 - length(matched_nonspiking)/(0.5*(length(non_spike_i) + length(non_spike_j)));
    end
end
clear i non_spike_i j non_spike_j matched_nonspiking

avg_nonspike_sim = mean(non_spiking_similarity(non_spiking_similarity ~= 1),'all');
std_nonspike_sim = std(non_spiking_similarity(non_spiking_similarity ~= 1),[],'all');

%Spike count similarity
average_spikes = zeros(1,clusters);
std_spikes = zeros(1,clusters);
for i = 1:clusters
    average_spikes(1,i) = mean(all_cluster_data(i).spike_count);
    std_spikes(1,i) = std(all_cluster_data(i).spike_count);
end
clear i
%plot of average spike counts
figure;
scatter([1:clusters],average_spikes)
hold on
errorbar([1:10],average_spikes,std_spikes,'vertical','LineStyle','none')
title('Average Number of Spikes per Simulation','FontSize',16)
xlabel('Initially Spiking Cluster Number','FontSize',16)
ylabel('Number of Spikes','FontSize',16)
xlim([0,clusters+1])
%calculate t-test of spike counts (only for spiking neurons
spike_count_h = zeros(clusters,clusters);
spike_count_p = zeros(clusters,clusters);
spike_sig_h = zeros(clusters,clusters);
for i = 1:clusters
    spiking_i = all_cluster_data(i).spike_count;
    for j = i+1:clusters
        spiking_j = all_cluster_data(j).spike_count;
        [p,h] = ranksum(spiking_i,spiking_j);
        spike_count_h(i,j) = h;
        spike_count_h(j,i) = h;
        spike_count_p(i,j) = p;
        spike_count_p(j,i) = p;
        spike_sig_h(i,j) = ansaribradley(spiking_i,spiking_j);
        spike_sig_h(j,i) = ansaribradley(spiking_i,spiking_j);
    end
end
clear i spiking_i j spiking_j h p
avg_num_from_diff_distributions = mean(spike_count_h,'all');

%Visualize which neurons spike for which initialization:
neuron_activity_per_initialization = zeros(clusters,n);
for i = 1:clusters
    spiking = setdiff(neuron_reorder,all_cluster_data(i).non_spiking);
    first_spikers = find(cluster_mat(i,:));
    neuron_activity_per_initialization(i,spiking) = neuron_activity_per_initialization(i,spiking) + 1;
    neuron_activity_per_initialization(i,first_spikers) = neuron_activity_per_initialization(i,first_spikers) + 1;
end
figure;
imagesc(neuron_activity_per_initialization)
colormap(jet)
colorbar('Ticks',[0,1,2],'TickLabels',{'Non-spiking','Spiking','Initially Spiking'})
title('Visualization of Neuronal Activity Per Simulation','FontSize',16)
ylabel('Simulation Number','FontSize',16)
xlabel('Neuron Number','FontSize',16)

%Neuron type similarity
excitatory_spike_sim = zeros(1,clusters);
inhibitory_spike_sim = zeros(1,clusters);
excitatory_nonspike_sim = zeros(1,clusters);
inhibitory_nonspike_sim = zeros(1,clusters);
for i = 1:clusters
    nonspike_clust = all_cluster_data(i).non_spiking;
    spike_clust = setdiff(neuron_reorder(:,i),nonspike_clust);
    E_clust_spike = intersect(spike_clust,E_indices);
    I_clust_spike = intersect(spike_clust,I_indices);
    E_clust_nonspike = intersect(nonspike_clust,E_indices);
    I_clust_nonspike = intersect(nonspike_clust,I_indices);
    excitatory_spike_sim(i) = length(E_clust_spike)/length(spike_clust);
    inhibitory_spike_sim(i) = length(I_clust_spike)/length(spike_clust);
    excitatory_nonspike_sim(i) = length(E_clust_nonspike)/length(nonspike_clust);
    inhibitory_nonspike_sim(i) = length(I_clust_nonspike)/length(nonspike_clust);
end
clear i nonspike_clust spike_clust E_clust_spike I_clust_spike ...
    E_clust_nonspike I_clust_nonspike

avg_percent_E_spike = mean(excitatory_spike_sim);
avg_percent_I_spike = mean(inhibitory_spike_sim);
avg_percent_E_nonspike = mean(excitatory_nonspike_sim);
avg_percent_I_nonspike = mean(inhibitory_nonspike_sim);
all_avg = [avg_percent_E_spike, avg_percent_I_spike, ...
    avg_percent_E_nonspike, avg_percent_I_nonspike];
std_percent_E_spike = std(excitatory_spike_sim);
std_percent_I_spike = std(inhibitory_spike_sim);
std_percent_E_nonspike = std(excitatory_nonspike_sim);
std_percent_I_nonspike = std(inhibitory_nonspike_sim);
all_std = [std_percent_E_spike, std_percent_I_spike, ...
    std_percent_E_nonspike, std_percent_I_nonspike];
%Scatter of average values w/ errorbars
figure;
scatter([1:4],all_avg)
hold on
errorbar([1:4],all_avg,all_std,'vertical','LineStyle','none')
title('Average Simulation Neuron Quality','FontSize',16)
xlabel('Condition','FontSize',16)
ylabel('Fraction of Simulation Population','FontSize',16)
xlim([0,5])
xticks([1,2,3,4])
xticklabels({'Excitatory Spiking Neurons', 'Inhibitory Spiking Neurons', ...
    'Excitatory Nonspiking Neurons', 'Inhibitory Nonspiking Neurons'})
%Plot of average values
figure;
plot(excitatory_spike_sim)
hold on
plot(inhibitory_spike_sim)
plot(excitatory_nonspike_sim)
plot(inhibitory_nonspike_sim)
legend({'Excitatory Spiking Fraction', 'Inhibitory Spiking Fraction', 'Excitatory Non-Spiking Fraction', 'Inhibitory Non-Spiking Fraction'})
title('Average Neuron Quality Per Simulation','FontSize',16)
xlabel('Simulation Number','FontSize',16)
ylabel('Fraction of Simulation Population','FontSize',16)
%Check that the populations significantly match the general population
%all h must be 0 for all to be Gaussian
h = adtest(excitatory_spike_sim);
h = adtest(inhibitory_spike_sim);
h = adtest(excitatory_nonspike_sim);
h = adtest(inhibitory_nonspike_sim);
%compare to a normal distribution with p_E or p_I mean and unknown sigma via ttest
rand_E = (mean(all_std))*randn([1,clusters])+p_E;
rand_I = (mean(all_std))*randn([1,clusters])+p_I;
ttest_E_spike = ttest(excitatory_spike_sim, rand_E);
ttest_I_spike = ttest(inhibitory_spike_sim, rand_I);
ttest_E_nonspike = ttest(excitatory_nonspike_sim, rand_E);
ttest_I_nonspike = ttest(inhibitory_nonspike_sim, rand_I);

%% Test conditions for regular spiking
%This code block will run through varying values of certain variables and
%store spiking data in a matrix. This will allow users to determine what
%conditions are necessary to create periodic spiking.

%Variables
%   We'll be looking at the synaptic conductance step size and the
%   probability of excitatory neurons in particular to determine a phase
%   plane of activity in this network
num_tests = 20;
del_G_syn_E_vec = linspace(1*10^(-9),20*10^(-9),num_tests);
del_G_syn_I_vec = linspace(1*10^(-9),20*10^(-9),num_tests);

avg_firing_rate = zeros(num_tests,num_tests); %rows are del_G_syn, columns are p_E
avg_total_spikes = zeros(num_tests,num_tests);
avg_ISI = zeros(num_tests,num_tests);

type = 'cluster'; %'neuron';
reset = 1;

t_start = 1;
t_end = t_steps;

for i = 1:length(del_G_syn_E_vec)
    display(strcat('Step ',string(i)))
    del_G_syn_E = del_G_syn_E_vec(i);
    for j = 1:length(del_G_syn_I_vec)
        del_G_syn_I = del_G_syn_I_vec(j);
        [V_m, G_sra, G_syn_I, G_syn_E, I_syn, conns] = lif_sra_calculator(n, ...
            clusters,cluster_mat,conns,conn_bin,V_m,V_reset,V_th,G_sra,del_G_sra,G_syn_I, ...
            G_syn_E, del_G_syn_E,del_G_syn_I,I_syn,E_syn_E,E_syn_I,E_K,E_L,G_L,C_m,t_start,t_end,dt, ...
            tau_sra,tau_syn_E,tau_syn_I,connectivity_gain,type,reset,I_indices,E_indices);
        spikes_V_m = V_m >= V_th;
        spike_count = sum(spikes_V_m,2);
        %calculate the average total number of spikes (of spiking neurons)
        avg_total_spikes(i,j) = mean(spike_count(spike_count ~= 0),'all');
        %calculate the average firing rate (of spiking neurons)
        avg_firing_rate(i,j) = mean(t_max*(spike_count(spike_count ~= 0)/t_steps),'all'); %in seconds
        %find ISIs
        [firing_x,firing_t] = find(spikes_V_m);
        comb_firing = [firing_x,firing_t];
        comb_firing_sort = sortrows(comb_firing);
        k = comb_firing_sort(1,1);
        ISIs = [];
        for l = 2:length(comb_firing_sort)
            if comb_firing_sort(l,1) == k
                ISIs(end+1) = comb_firing_sort(l,2) - comb_firing_sort(l-1,2); %#ok<SAGROW>
            else
                k = comb_firing_sort(l,1);
            end    
        end
        ISIs = ISIs*dt; %convert to seconds
        ISIs(isnan(ISIs)) = [];
        mean_ISI = mean(ISIs);
        if isnan(mean_ISI) %empty
            avg_ISI(i,j) = inf;
        else
            avg_ISI(i,j) = mean(ISIs);
        end
    end
end
%Visualize
%Average Firing Rate
figure;
imagesc(avg_firing_rate)
colorbar
xticks([1:num_tests])
xticklabels(del_G_syn_I_vec)
xlabel('\Delta G_{synI}','FontSize',20)
yticks([1:num_tests])
yticklabels(del_G_syn_E_vec)
ylabel('\Delta G_{synE}','FontSize',20)
title('Average Firing Rate (s)','FontSize',20)
%Average Total Spikes
figure;
imagesc(avg_total_spikes)
colorbar
xticks([1:num_tests])
xticklabels(del_G_syn_I_vec)
xlabel('\delta G_{syn I}','FontSize',20)
yticks([1:num_tests])
yticklabels(del_G_syn_E_vec)
ylabel('\delta G_{syn E}','FontSize',20)
title('Average Total Spikes','FontSize',20)
%Average ISI (in seconds)
figure;
imagesc(avg_ISI)
colorbar
xticks([1:num_tests])
xticklabels(del_G_syn_I_vec)
xlabel('\delta G_{syn I}','FontSize',20)
yticks([1:num_tests])
yticklabels(del_G_syn_E_vec)
ylabel('\delta G_{syn E}','FontSize',20)
title('Average ISI (s)','FontSize',20)

filename = '/Users/hannahgermaine/Documents/PhD/Rotations/Paul Miller/Code/lif_sra_replay_attractor/test_results/del_syn_tests/';
save(strcat(filename,'avg_firing_rate','_n',string(n),'_cp',strrep(string(conn_prob),'.','_')),'avg_firing_rate')
save(strcat(filename,'avg_total_spikes','_n',string(n),'_cp',strrep(string(conn_prob),'.','_')),'avg_total_spikes')
save(strcat(filename,'avg_ISI','_n',string(n),'_cp',strrep(string(conn_prob),'.','_')),'avg_ISI')

%% Visualize results of above tests

filename = '/Users/hannahgermaine/Documents/PhD/Rotations/Paul Miller/Code/lif_sra_replay_attractor/test_results/del_syn_tests/';
files = dir(filename);
files = files([files.isdir] == 0);
for i = 1:length(files)
    dataname = files(i).name;
    data = load(strcat(filename,dataname));
    name = fieldnames(data);
    figure(i);
    imagesc(data.(name{1}))
    title(dataname)
    colorbar
end    

%% Create a Network Visualization

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