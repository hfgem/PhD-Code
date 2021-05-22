function [V_m, G_sra, G_syn_I, G_syn_E, I_syn, conns] = lif_sra_calculator(n, ...
    clusters,cluster_mat,conns,V_m,V_reset,V_th,G_sra,del_G_sra,...
    G_syn_I, G_syn_E, del_G_syn_E,del_G_syn_I,I_syn,E_syn_E,E_syn_I,E_K,E_L,...
    G_L,C_m,t_start,t_end,dt,tau_sra,tau_syn_E,tau_syn_I,connectivity_gain,...
    type,reset,I_indices,E_indices,I_theta_in,seed)

    rng(seed)
    %Assumptions:
    %   1. LIF model with spike rate adaptation and synaptic transmission
    %   2. SRA has decay
    %   3. Synaptic transmission has decay
    %   4. Strength of connections from presynaptic neurons is equal for 
    %       all postsynaptic neurons
    %   5. There can be an externally applied current input through 
    %       I_theta_in
    %   6. Excitatory and inhibitory connections have different reversal
    %       potentials in the postsynaptic neuron represented in a matrix
    %       E_syn with rows representing the presynaptic neuron, and
    %       columns the postsynaptic neuron
    %   6. Excitatory and inhibitory currents have different synaptic
    %       conductance steps and decay rates
    %   7. There is no noise added at each timestep to the final membrane
    %       potential calculation.
    %   8. delta G_syn is proportional to the connection strength


    %Determine whether the first spikers are individual neurons / a cluster
    %and pick which cluster or neurons are initially spiking
    if strcmp(type,'cluster')
        i = randi(clusters);
        display(i) %if you want to see which cluster was chosen, uncomment
        neur_start = find(cluster_mat(i,:)); %one cluster starts firing
        clear i
    else
        neur_start = randi(n,[1,clusters]); %cluster # of neurons start firing
    end
    
    %Run any reset settings
    if (t_start == 1) || (reset == 1) %reset firing rate to 0 if first or reset
        V_m(:,t_start) = V_reset; %set all neurons to baseline reset membrane potential with added noise
        G_sra(:,t_start) = 0;
        G_syn_E(:,t_start) = 0;
        G_syn_I(:,t_start) = 0;
    end
    
    %Set the membrane potential of spiking neurons to threshold
%     V_m(neur_start,t_start) = V_th;
    %AND/OR (experimental) set the input current to only apply to the
    %'initiator' neurons
    I_theta = ones(size(V_m));
%     I_input(neur_start,:) = 1;
    I_theta = I_theta_in.*I_theta;
    
    %Run through each timestep and calculate
    for t = t_start:t_end
        %check for spiking neurons and postsynaptic and separate into E and I
        spikers = find(V_m(:,t) >= V_th);
        spikers_I = spikers(ismember(spikers,I_indices)); %indices of inhibitory presynaptic neurons
        spikers_E = spikers(ismember(spikers,E_indices)); %indices of excitatory presynaptic neurons
        %______________________________________
        %Adjust parameters dependent on spiking
        G_sra(spikers,t) = G_sra(spikers,t) + del_G_sra; %set SRA conductance values
        %Synaptic conductance is stepped for postsynaptic neurons
        %   dependent on the number of presynaptic connections, and the
        %   current will depend on the presynaptic neuron type (E_syn_I and E_syn_E)
        incoming_conn_E = sum(conns(spikers_E,:),1)';
        incoming_conn_I = sum(conns(spikers_I,:),1)';
        G_syn_E(:,t) = G_syn_E(:,t) + del_G_syn_E*incoming_conn_E; %excitatory conductance
        G_syn_I(:,t) = G_syn_I(:,t) + del_G_syn_I*incoming_conn_I; %inhibitory conductance
        I_syn(:,t) = G_syn_E(:,t).*(E_syn_E - V_m(:,t)) + G_syn_I(:,t).*(E_syn_I - V_m(:,t));
        I_app = I_syn(:,t) + I_theta(:,t);
        %______________________________________
        %Calculate membrane potential using integration method
        V_ss = (I_app + G_L*E_L + G_sra(:,t)*E_K)./(G_L + G_sra(:,t)); %"steady state" calculation
        taueff = C_m./(G_L + G_sra(:,t)); %timescale for change in the membrane potential
        exp_coeff = (G_L*(E_L - V_m(:,t)) + G_sra(:,t).*(E_K - V_m(:,t)) + I_app)./(G_L + G_sra(:,t));
        V_m(:,t+1) = V_ss - exp_coeff.*exp(-dt./taueff) + randn([n,1])*(10^(-4));
        V_m(spikers,t+1) = V_reset; %update those that just spiked to reset
        %______________________________________
        %Update next step conductances
        G_sra(:,t+1) = G_sra(:,t)*exp(-dt/tau_sra); %Spike rate adaptation conductance
        %Synaptic conductance updated for each postsynaptic neuron by
        %incoming connection type
        G_syn_E(:,t+1) = G_syn_E(:,t).*exp(-dt/tau_syn_E); %excitatory conductance update
        G_syn_I(:,t+1) = G_syn_I(:,t).*exp(-dt/tau_syn_I); %inhibitory conductance update
        %______________________________________
        %Update connection strengths
        conns(spikers,:) = connectivity_gain*conns(spikers,:); %enhance connections of those neurons that just fired
        conns(~spikers,:) = (1/connectivity_gain)*conns(~spikers,:); %decrease connections of those neurons that did not fire
    end
end