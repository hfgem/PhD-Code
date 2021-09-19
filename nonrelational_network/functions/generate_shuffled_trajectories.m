function [shuffled_spike_sequences] = generate_shuffled_trajectories(n,...
    shuffle_n, network_spike_sequences, viable_inits)
    %_________
    %ABOUT: This function generates shuffled firing sequences based on the
    %statistics of real trajectories from network simulations.
    %
    %INPUTS:
    %   n = number of neurons in the simulation
    %   shuffle_n = number of shuffled trajectories to generate
    %   network_spike_sequences = a structure file with real spike
    %       sequences generated by the network. Fields contained are as
    %       follows:
    %           1. events = vector of event start and stop times
    %           2. spike_order = order of spiking neurons, excluding
    %               nonspiking neurons
    %           3. spike_ranks = vectors for each event of each neuron's
    %               rank in the spike order, with '0's for neurons that do
    %               not spike
    %           4. nonspiking_neurons = binary vectors for each event of
    %               which neurons did not spike (1) and which did (0).
    %   viable_inits = vector of indices of viable initializations in the 
    %       simulation data
    %
    %OUTPUTS:
    %   shuffled_spike_sequences = a structure file with generated spike
    %       sequences following the statistics of the network spike
    %       sequences. Fields contained are as follows:
    %           1. spike_order = order of spiking neurons, excluding
    %               nonspiking neurons
    %           2. spike_ranks = vectors for each event of each neuron's
    %               rank in the spike order, with '0's for neurons that do
    %               not spike
    %           3. nonspiking_neurons = binary vectors for each event of
    %               which neurons did not spike (1) and which did (0).
    %_________
    
    %First gather generated spike sequences statistics
    spiking_nums = [];
    for i = viable_inits
        if ~isempty(network_spike_sequences(i).spike_ranks) %second check
            sequences = network_spike_sequences(i).spike_order(1);
            sequence_names = fieldnames(sequences);
            [num_seq,~] = size(sequence_names);
            for j = 1:num_seq
                number_neurons = length(network_spike_sequences(i).spike_order.(sequence_names{j}));
                if number_neurons >= 0.25*n %at least 1/4 of the neurons participate in the sequence
                    spiking_nums = [spiking_nums; number_neurons]; %#ok<AGROW>
                end
            end  
        end
    end   
    avg_num_spiking = mean(spiking_nums);
    std_num_spiking = std(spiking_nums);
    
    %Next generate shuffled data based on spike statistics
    shuffled_spike_sequences = struct;
    for i = 1:shuffle_n
        %generate spike sequence randomly
        num_gen = round(avg_num_spiking + std_num_spiking*randn());
        if num_gen > n
            num_gen = n;
        end
        shuffle_seq = randperm(n,num_gen);
        shuffled_spike_sequences(i).spike_order.sequence_1 = shuffle_seq;
        %store ranks for each neuron
        ranks_vec = zeros(1,n);
        for j = 1:length(shuffle_seq)
            n_ind = shuffle_seq(j);
            ranks_vec(1,n_ind) = j;
        end
        shuffled_spike_sequences(i).spike_ranks.sequence_1 = ranks_vec;
        %store nonspiking neurons
        nonspiking_neurons = isnan(ranks_vec./ranks_vec);
        shuffled_spike_sequences(i).nonspiking_neurons.sequence_1 = nonspiking_neurons;
    end    
    
end