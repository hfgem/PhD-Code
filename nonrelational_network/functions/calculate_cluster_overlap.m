function [cluster_prog] = calculate_cluster_overlap(n,clusters,init,V_m_struct, ...
    bin_frac, events)
    %_________
    %ABOUT: This function calculates the cluster overlap of neurons spiking
    %in each event in each initialization for a particular network setup.
    %
    %INPUTS:
    %   n = number of neurons in the simulation
    %   clusters = number of clusters in the simulation
    %   init = number of initializations in the simulation
    %   V_m_struct = a structure variable with both the membrane potential
    %       array (V_m) and the binary spike array (spikes_V_m) for each
    %       initialization.
    %   bin_frac = what percentage of the total event length should be used
    %       in a single bin size.
    %   event = a [1 x init] cell with the event start and end times for
    %       all events in each initialization.
    %
    %OUTPUTS:
    %   cluster_prog = a structure file with cluster overlap statistics for
    %       each bin in each event. The statistics include the bin size,
    %       the number of clusters represented in each bin, the number of
    %       neurons represented in each bin, and the cluster overlap
    %       fraction.
    %_________
    
    cluster_prog = struct; %create structure to save all data
    for i = 1:init
        %Pull spike data
        spikes_V_m = V_m_struct(i).spikes_V_m;
        [event_num, ~] = size(events{i});
        event_times = events{i};
        for j = 1:event_num %for each event calculate cluster overlap
            event_length = event_times(j,2) - event_times(j,1);
            event_V_m = spikes_V_m(:,event_times(j,1):event_times(j,2));
            [~,spikes_t] = find(event_V_m);
            total_spikes = length(unique(spikes_t));
            cluster_prog(i).event(j).total_spikes = total_spikes;
            if total_spikes > 1 %First make sure we have multiple spikes to work with
                bin_size = ceil(event_length*bin_frac);
                %Run through each dataset and calculate overlap fractions
                cluster_number = [];
                neuron_number = [];
                cluster_overlap = [];
                for k = 1:event_length - bin_size %slide bin along all event time
                    spike_cut = event_V_m(:,j:j+bin_size);
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
                clear k spike_cut n_spiking c_spiking n_c avg_overlap
                %Store relevant information
                cluster_prog(i).event(j).bin_size = bin_size;
                cluster_prog(i).event(j).cluster_number = cluster_number;
                cluster_prog(i).event(j).neuron_number = neuron_number;
                cluster_prog(i).event(j).cluster_overlap = cluster_overlap;
                clear bin_size cluster_number neuron_number cluster_overlap
            end 
            clear event_length event_V_m spikes_t total_spikes
        end
        clear spikes_V_m event_num event_times j
    end
    clear i
end