function [neuron_ranks, rho_unique, ranks, rho_nonunique, ranks_mod] = calculate_trajectory_similarity(n, ...
    inits, network_spike_sequences, network_spike_ranks, non_spiking_neurons)
    %_________
    %INPUTS:
    %   n = number of neurons in the simulation
    %   inits = number of initializations in the simulation
    %   network_spike_sequences = a [inits x n] size matrix of the order
    %           in which neurons spiked
    %   network_spike_ranks = a [inits x n] size matrix of the rank of
    %           each neuron, with nonspiking neurons marked as 0 rank
    %   non_spiking_neurons = a [inits x n] size binary matrix of which
    %           neurons are non-spiking in each particular initialization
    %
    %OUTPUTS:
    %   neuron_ranks = a modified version of network_spike_ranks that ranks
    %           all nonspiking neurons at the end in order of index.
    %   rho_unique = a [inits x inits] matrix of Spearman's rank correlation
    %           rhos for spike ranks with the nonspiking neurons ranked at  
    %           the end of the list in order of index - using Matlab's 
    %           function
    %   ranks = a [inits x inits] size matrix of manually calculated 
    %           correlation ranks between pairs of initializations
    %   rho_nonunique = a [inits x inits] size matrix of Spearman's
    %           correlation rhos for spike ranks excluding nonspiking
    %           neurons using Matlab's function
    %   ranks_mod = a [inits x inits] size matrix of manually calculated
    %           correlation ranks excluding nonspiking neurons. Neurons
    %           that are nonspiking in one sequence but not another are
    %           kept, and neurons that are nonspiking in both are excluded.
    %_________
    
    %======NEURON RANKS - INCLUDING NONSPIKING======

    %calculate ranks (including nonspiking)
    neuron_ranks = zeros(n,inits);
    for i = 1:inits
        for j = 1:n
            ind = find(network_spike_sequences(i,:) == j);
            neuron_ranks(j,i) = ind;
        end
    end
    clear i j ind

    %Using Spearman's rank correlation on the unique integer values we get
    rho_unique = corr(neuron_ranks,'Type','Spearman');

    %If all n ranks are distinct integers, can use the formula: (wikipedia)
    %r_s = 1 - (6*sum(d_i^2))/(n*(n^2 - 1))
    %d_i = rg(X_i) - rg(Y_i)
    ranks = zeros(inits,inits);
    for i = 1:inits
        X_i = neuron_ranks(:,i);
        for j = 1:inits
            Y_i = neuron_ranks(:,j);
            d = zeros(1,n);
            for k = 1:n
                ix_i = find(X_i == k);
                iy_i = find(Y_i == k);
                d(1,k) = ix_i - iy_i;
            end
            ranks(i,j) = 1 - (6*sum(d.^2))/(n*(n^2 - 1));
            if i == j
                ranks(i,j) = 1;
            end
        end
    end
    clear i X_i Y_i d j ix_i iy_i k
    
    %======NEURON RANKS - NOT INCLUDING NONSPIKING======

    %Rank calculation with all nonspiking neurons ranked as NaN
    network_spike_ranks_nan = network_spike_ranks;
    network_spike_ranks_nan(network_spike_ranks_nan == 0) = NaN;
    rho_nonunique = corr(network_spike_ranks_nan','Type','Spearman','Rows','pairwise');

    %Rank calculation modified to exclude neurons that are nonspiking in
    %both datasets, but still including nonspiking neurons that are in one
    %ordered at the end of the dataset
    ranks_mod = zeros(inits,inits);
    for i = 1:inits
        X_i = neuron_ranks(:,i);
        X_nonspike = find(non_spiking_neurons(i,:));
        for j = 1:inits
            Y_i = neuron_ranks(:,j);
            Y_nonspike = find(non_spiking_neurons(j,:));
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
                if ~isempty(ix_i) && ~isempty(iy_i)
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

end