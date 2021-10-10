%Testing how STDP applied to one sequence affects other sequences

%______ABOUT______:
%This code tests how repeating one sequence through re-initialization
%affects other sequences. The code also makes use of previous network
%structures and parameters that generated successful sequences (from
%lif_network_postrotation.m) to load a constant network structure for
%testing.
%__________________

%% Initialization

%Select folder of good network structure and parameters to use in tests
load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select network.m Load Folder'); %Have user input where they'd like network structure to be loaded from

%_________________________________
%___Load Independent Parameters___
%_________________________________
load(strcat(load_path,'/network.mat'))
slashes = find(load_path == '/');
param_path = load_path(1:slashes(end));
load(strcat(param_path,'/parameters.mat'))

%TO DEPRECATE LATER: In the case that the parameters file being used 
%does not contain V_m_noise:
try
    parameters.V_m_noise;
catch
    parameters.('V_m_noise') = 10^(-4);
end

%TO DEPRECATE LATER: In the case that the parameters file being used 
%does not contain tau_stdp:
try
    parameters.tau_stdp;
catch
    parameters.('tau_stdp') = 5*10^(-3);
end

%____________________________________
%___Calculate Dependent Parameters___
%____________________________________
cluster_n = min(parameters.n*2/parameters.clusters,parameters.n); %number of neurons in a cluster (for small n round(n/3), for large n round(n/5)) 
parameters.('cluster_n') = cluster_n;

%Interaction constants
t_steps = parameters.t_max/parameters.dt; %number of timesteps in simulation
E_syn_E = parameters.V_syn_E*ones(parameters.n,1); %vector of the synaptic reversal potential for excitatory connections
E_syn_I = parameters.V_syn_I*ones(parameters.n,1); %vector of the synaptic reversal potential for inhibitory connections
IES = ceil(parameters.IEI/parameters.dt); %inter-event-steps = the number of steps to elapse between spikes
%save for easy calculations
parameters.('t_steps') = t_steps;
parameters.('E_syn_E') = E_syn_E;
parameters.('E_syn_I') = E_syn_I;
parameters.('IES') = IES;

%How many tests of different initializations to run
if strcmp(parameters.type,'cluster')
    test_val_max = parameters.clusters; %Ensures every cluster is initialized
else
    test_val_max = 10; %This value can be modified as you see fit
end
%save for easy calculations
parameters.('test_val_max') = test_val_max;

%Adding an input current to all cells (one of the options must be uncommented)
x_in = [0:parameters.dt:parameters.t_max];
% %Rhythmic current input: (uncomment if desired)
% I_coeff = 0; %5.1*10^(-10); %set to 0 for no input current
% I_Hz = 1; %frequency of input - 1 Hz or 60 bpm, 4-8 Hz for theta rhythm
% I_in = I_coeff*(0.5 + 0.5*sin(I_Hz*2*pi*x_in)); %Approximately theta wave input current
%Noisy input: (uncomment if desired)
I_in = parameters.I_coeff*randn(parameters.n,parameters.t_steps+1)*parameters.I_scale; %Generally use -0.5-0.5 nA stimulus
%save for easy calculations
parameters.('x_in') = x_in;
parameters.('I_in') = I_in;

%Calculate connection probabilites
npairs = parameters.n*(parameters.n-1); %total number of possible neuron connections
nclusterpairs = parameters.cluster_n*(parameters.cluster_n - 1)*parameters.clusters; %total number of possible intra-cluster connections
cluster_prob = min(parameters.conn_prob*npairs/nclusterpairs,1); %0.2041; %intra-cluster connection probability
p_I = 1 - parameters.p_E; %probability of an inhibitory neuron
n_I = round(p_I*parameters.n); %number of inhibitory neurons
%save for easy calculations
parameters.('npairs') = npairs;
parameters.('nclusterpairs') = nclusterpairs;
parameters.('cluster_prob') = cluster_prob;
parameters.('p_I') = p_I;
parameters.('n_I') = n_I;

%____________________________________
%___Assign STDP Rules Desired________
%____________________________________
%Here you can modify the STDP rules previously saved
parameters.type = 'neuron'; %Each sequence initialization will be through setting neurons to threshold
parameters.tau_stdp = 1*10^(-3); %STDP time constant (s)
parameters.connectivity_gain = 0.002; %amount to increase or decrease connectivity by with each STDP rule (more at the range of 0.002-0.005)
%Here you set the simulation rules
num_repeat_tests = 10; %how many repeat values to test
max_repeats = 100; %maximum number of repeats to test
num_repeats = 0:num_repeat_tests:max_repeats; %vector of repeats to test
num_repeats(1) = 1; %update first value to be 1 repeat rather than 0

%% Test STDP's Effect on Other Sequences

%Run learning on one of the initializations the number of times
%specified in num_repeats and save the final sequence of the
%initialization used, and then run the initializations for other sequences
%and save them

%Select simulation save folder
save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder'); %Have user input where they'd like network structure to be loaded from

%Create a structure file to store all the sequence results
STDP_sequences = struct;
STDP_V_m = struct;

for i = 1:num_repeat_tests+1
    rng(1) %to ensure the same initialization is selected each time
    init_seed = randi(parameters.test_val_max,[1, parameters.test_val_max]);
    
     %SET UP NETWORK
    [cluster_mat, conns] = create_clusters(parameters, i, 1);
    conns_copy = conns; %just a copy of the connections to maintain for reset runs if there's "plasticity"
    %Randomize excitatory and inhibitory connection strengths based on selected
    %probability.
    all_indices = [1:parameters.n];
    I_indices = datasample(all_indices,parameters.n_I,'Replace',false); %indices of inhibitory neurons
    E_indices = find(~ismember(all_indices,I_indices)); %indices of excitatory neurons
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
    
    %Set up storage
    STDP_sequences(i).num_repeats = num_repeats(i);
    STDP_sequences(i).init_seed = init_seed;
    STDP_sequences(i).spike_order = zeros(parameters.n,parameters.test_val_max);
    STDP_sequences(i).nonspiking = zeros(parameters.n,parameters.test_val_max);
    
    %Run for selected number of repeats
    %init_seed(1) will be the one that is repeated
    for j = 1:num_repeats(i)
        %Create Storage Variables
        I_syn = zeros(parameters.n,parameters.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
        %synaptic conductances for each neuron at each timestep
        G_syn_I = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic inhibitory (S)
        G_syn_E = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic excitatory (S)
        V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
        G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep
        
        %Run simulation
        [V_m, ~, ~, ~, ~, conns_update] = lif_sra_calculator_postrotation(...
            parameters, init_seed(1), network_copy, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
        %update network file with new connectivity matrix from STDP
        network_copy.conns = conns_update;
    end
    STDP_sequences(i).conns = network_copy.conns; %Save the final connectivity matrix
    STDP_V_m(i).(strcat('V_m_',string(1))) = V_m;
    
    %Store sequence from init_seed(1)
    spikes_V_m = V_m >= parameters.V_th;
    [spikes_x,~] = find(spikes_V_m);
    spike_order = unique(spikes_x,'stable');
    nonspiking = sum(spikes_V_m,2) == 0;
    STDP_sequences(i).spike_order(:,1) = [spike_order; find(nonspiking)];
    STDP_sequences(i).nonspiking(:,1) = nonspiking;
    
    %Clear unnecessary variables
    clear j I_syn G_syn_I G_syn_E V_m G_sra conns_update spikes_V_m ...
        spikes_x spike_order nonspiking
    
    %Run other sequence initializations
    for j = 2:parameters.test_val_max
        %Create Storage Variables
        I_syn = zeros(parameters.n,parameters.t_steps+1); %synaptic current emitted by each neuron at each timestep (A)
        %synaptic conductances for each neuron at each timestep
        G_syn_I = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic inhibitory (S)
        G_syn_E = zeros(parameters.n,parameters.t_steps+1); %conductance for presynaptic excitatory (S)
        V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3)); %set all neurons to baseline reset membrane potential with added noise
        G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep
        
        %Run simulation
        [V_m, ~, ~, ~, ~, ~] = lif_sra_calculator_postrotation(...
            parameters, init_seed(j), network_copy, I_syn, G_syn_I, G_syn_E, V_m, G_sra);
        STDP_V_m(i).(strcat('V_m_',string(j))) = V_m;
        
        %Store sequence from this initialization
        spikes_V_m = V_m >= parameters.V_th;
        [spikes_x,~] = find(spikes_V_m);
        spike_order = unique(spikes_x,'stable');
        nonspiking = sum(spikes_V_m,2) == 0;
        STDP_sequences(i).spike_order(:,j) = [spike_order; find(nonspiking)];
        STDP_sequences(i).nonspiking(:,j) = nonspiking;
        
        %Clear unnecessary variables
        clear j I_syn G_syn_I G_syn_E V_m G_sra conns_update spikes_V_m ...
            spikes_x spike_order nonspiking
    end 
end 

%SAVE DATA
save(strcat(save_path,'/STDP_sequences.mat'),'STDP_sequences','-v7.3')
save(strcat(save_path,'/STDP_V_m.mat'),'STDP_V_m','-v7.3')
save(strcat(save_path,'/network.mat'),'network_copy','-v7.3') %contains the updated connectivity matrix

disp('STDP Tests Done')

%% Analyze the Results of STDP on One Sequence
%Load data
%Select folder of STDP datasets
load_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder');
load(strcat(load_path,'/STDP_sequences.mat'))
load(strcat(load_path,'/STDP_V_m.mat'))
slashes = find(load_path == '/');
param_path = load_path(1:slashes(end));
load(strcat(param_path,'/parameters.mat'))

%_____________
% Visualize how the sequence changes over the course of STDP learning
%_____________
%Pull /data
num_repeats = [STDP_sequences.num_repeats];
num_to_visualize = length(num_repeats);
same_sequences = zeros(parameters.n,num_to_visualize);
same_non_spiking = zeros(parameters.n,num_to_visualize);
for i = 1:num_to_visualize
    same_sequences(:,i) = STDP_sequences(i).spike_order(:,1);
    same_non_spiking(:,i) = STDP_sequences(i).nonspiking(:,1);
end
%Set up plot
subplot_x = floor(sqrt(num_to_visualize));
subplot_y = ceil(sqrt(num_to_visualize));
f = figure;
axes = [];
%To order all sequences by a particular firing sequence, uncomment the
%following and select which sequence to order by by changing the y-index
% STDP_order = same_sequences(:,1);
for i = 1:num_to_visualize
    ax = subplot(subplot_x,subplot_y,i);
    V_m_spikes = STDP_V_m(i).('V_m_1') >= parameters.V_th;
    STDP_order = same_sequences(:,i); %Comment out if ordering by line 228
    V_m_spikes_reordered = V_m_spikes(STDP_order,:);
    spike_times = find(sum(V_m_spikes_reordered,1));
    large_spike_intervals = find((spike_times(2:end) - spike_times(1:end-1))>500);
    last_spike_time = spike_times(large_spike_intervals(1))+100;
    imagesc(V_m_spikes_reordered(:,1:last_spike_time))
    colormap(flip(gray))
    xlabel('Time (s)')%,'FontSize',16)
    ylabel('Reordered Neuron Number')%,'FontSize',16)
    title(strcat('Repeats = ',string(num_repeats(i))))
    axes = [axes, ax]; %#ok<AGROW>
end
linkaxes(axes)
sgtitle('Sequence as a Function of Learning')
savefig(f,strcat(load_path,'/sequences_after_learning.fig'))
saveas(f,strcat(load_path,'/sequences_after_learning.jpg'))
    

%_____________
% Calculate the sequence correlations for different amounts of learning
% This calculation excludes nonspiking neurons
%_____________
ranks = zeros(num_to_visualize,num_to_visualize);
for i = 1:num_to_visualize
    X_i = same_sequences(:,i);
    X_nonspike = find(same_non_spiking(:,i));
    for j = 1:num_to_visualize
        Y_i = same_sequences(:,j);
        Y_nonspike = find(same_non_spiking(:,j));
        nonspiking = unique([X_nonspike',Y_nonspike']);
        %Only look at neurons spiking in both
        X_i = setdiff(X_i, nonspiking', 'stable');
        Y_i = setdiff(Y_i, nonspiking', 'stable');
        new_n = (length(X_i)+length(Y_i))/2;
        d = zeros(1,parameters.n);
        for k = 1:parameters.n
            %find the ranks of each neuron
            ix_i = find(X_i == k);
            iy_i = find(Y_i == k);
            if ~isempty(ix_i) && ~isempty(iy_i) %if both are spiking
                d(1,k) = ix_i - iy_i;
            end
        end
        ranks(i,j) = 1 - (6*sum(d.^2))/(new_n*(new_n^2 - 1));
        ranks(j,i) = ranks(i,j);
        if i == j
            ranks(i,j) = 1;
        end
    end
end    
clear i X_i X_nonspike j Y_i Y_nonspike nonspiking new_n d k ix_i iy_i

ranks_vec = nonzeros(triu(ranks,1));

f1 = figure;
subplot(1,2,1)
histogram(ranks_vec,num_to_visualize)
title('Histogram of Unique Rank Correlation Values')
subplot(1,2,2)
imagesc(ranks)
colorbar()
xticks([1:num_to_visualize])
yticks([1:num_to_visualize])
xt = get(gca,'XTick');
yt = get(gca,'YTick');
set(gca, 'XTick',xt, 'XTickLabel',num_repeats)
set(gca, 'YTick',yt, 'YTickLabel',num_repeats)
title('Visualization of Rank Correlation Pairs')
savefig(f1,strcat(load_path,'/sequences_after_learning_correlation.fig'))
saveas(f1,strcat(load_path,'/sequences_after_learning_correlation.jpg'))

% Visualize histograms of correlation values for different learning amounts


% Visualize a line plot with std of the average correlation values for
% different learning amounts


