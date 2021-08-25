function [V_m, G_sra, G_syn_I, G_syn_E, I_syn, conns] = lif_sra_calculator_postrotation(...
    parameters, seed, cluster_mat, conns, I_indices, E_indices, I_syn, ...
    G_syn_I, G_syn_E, V_m, G_sra)
    %_________
    %ABOUT: This function uses the leaky integrate-and-fire model of 
    %neuronal firing to calculate the trajectory of membrane potentials,
    %currents, etc... that take place in a particular network with a 
    %particular set of parameters and initialization.
    %
    %INPUTS:
    %   parameters = a structure that contains the following:
    %       n = Number of neurons in the network
    %       V_reset = The reset membrane potential (V)
    %       V_th = The threshold membrane potential (V)
    %       
    %       del_G_sra = spike rate adaptation conductance step following spike 
    %               ranges from 1-200 *10^(-9) (S)
    %       del_G_syn_E = Synaptic conductance step for excitatory neurons 
    %               following spike (S)
    %       del_G_syn_I = Synaptic conductance step for inhibitory neurons 
    %               following spike (S)
    %       E_syn_E = An [n x 1] vector of the synaptic reversal potential for
    %               excitatory connections (V)
    %       E_syn_I = An [n x 1] vector of the synaptic reversal potential for
    %               inhibitory connections (V)
    %       E_K = Potassium reversal potential (V)
    %       E_L = Leak reversal potential (V)
    %       G_L = Leak conductance (S) - 10-30 nS range
    %       C_m = Total membrane capacitance (F)
    %       dt = Timestep (s)
    %       tau_sra = Spike rate adaptation time constant (s)
    %       tau_syn_E = AMPA/NMDA synaptic decay time constant (s)
    %       tau_syn_I = GABA synaptic decay time constant (s)
    %       connectivity_gain = Amount to increase or decrease connectivity by 
    %               with each spike (more at the range of 1.002-1.005) -
    %               keep at 1 to ensure no connectivity change
    %       I_in = An [n x t_steps + 1] matrix with input current values
    %       t_steps = The number of timesteps in the simulation
    %       type = Determines how spiking is initiated. Either:
    %           1. type = 'cluster', which sets a cluster of neurons to
    %           threshold at step 1
    %           2. type = 'neuron', which sets a fraction of neurons to
    %           threshold at step 1, where the fraction is set on line 
    %           3. type = 'current', which means spiking depends entirely on
    %           I_in and noise for initialization
    %   seed = A random number generator seed which:
    %       1. when type = 'cluster' sets which cluster is to be used for
    %           the initialization of spiking
    %       2. when type = 'neuron' sets the random number generator seed
    %           which affects the random selection of neurons used in the
    %           initialization of spiking as well as the random noise added
    %           to the membrane potential
    %       3. when type = 'current' sets the random number generator seed 
    %           which affects only the random noise added to the membrane 
    %           potential 
    %   cluster_mat = A binary [clusters x n] matrix of which neurons are
    %               in which cluster
    %   conns = An [n x n] matrix of which neurons are connected to each
    %               other, with values greater than 1 implying stronger
    %               connectivity
    %   I_indices = Vector of indices of inhibitory neurons
    %   E_indices = Vector of indices of excitatory neurons
    %   I_syn = An [n x t_steps+1] matrix of synaptic current emitted by 
    %               each neuron at each timestep (A)
    %   G_syn_I = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               inhibitory (S)
    %   G_syn_E = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               excitatory (S)
    %   V_m = An [n x t_steps+1] matrix of membrane potential for each 
    %               neuron at each timestep
    %   G_sra = An [n x t_steps+1] matrix with refractory conductance for
    %               each neuron at each timestep (S)
    %
    %OUTPUTS:
    %   V_m = An [n x t_steps+1] matrix of membrane potential for each 
    %               neuron at each timestep
    %   G_sra = An [n x t_steps+1] matrix with refractory conductance for
    %               each neuron at each timestep (S)
    %   G_syn_I = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               inhibitory (S)
    %   G_syn_E = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               excitatory (S)
    %   I_syn = An [n x t_steps+1] matrix of synaptic current emitted by 
    %               each neuron at each timestep (A)
    %   conns = The updated conns matrix (if connectivity_gain != 1)
    %
    %ASSUMPTIONS:
    %   1. LIF model with spike rate adaptation and synaptic transmission
    %   2. SRA has decay
    %   3. Synaptic transmission has decay
    %   4. There can be an externally applied current input through 
    %       I_in
    %   5. Excitatory and inhibitory connections have different reversal
    %       potentials in the postsynaptic neuron represented in vectors
    %       E_syn_E and E_syn_I
    %   6. Excitatory and inhibitory currents have different synaptic
    %       conductance steps and decay rates
    %_________

    %Pick which cluster is initially spiking
    if strcmp(parameters.type,'cluster') %cluster set to threshold at timestep 1
        rng(1)
        neur_start = find(cluster_mat(seed,:)); %one cluster starts firing
        %Set the membrane potential of spiking neurons to threshold
        V_m(neur_start,1) = parameters.V_th; %#ok<FNDSB>
    elseif strcmp(parameters.type,'neuron')
        rng(seed)
        neur_start = rand(parameters.n,1) <= 0.05;
        %Set the membrane potential of spiking neurons to threshold
        V_m(neur_start,1) = parameters.V_th; 
    elseif strcmp(parameters.type,'current')
        rng(seed)
    end
    I_theta = parameters.I_in;
    
    %Run through each timestep and calculate
    for t = 1:parameters.t_steps
        %check for spiking neurons and postsynaptic and separate into E and I
        spikers = find(V_m(:,t) >= parameters.V_th);
        spikers_I = spikers(ismember(spikers,I_indices)); %indices of inhibitory presynaptic neurons
        spikers_E = spikers(ismember(spikers,E_indices)); %indices of excitatory presynaptic neurons
        %______________________________________
        %Adjust parameters dependent on spiking
        G_sra(spikers,t) = G_sra(spikers,t) + parameters.del_G_sra; %set SRA conductance values
        %Synaptic conductance is stepped for postsynaptic neurons
        %   dependent on the number of presynaptic connections, and the
        %   current will depend on the presynaptic neuron type (E_syn_I and E_syn_E)
        incoming_conn_E = sum(conns(spikers_E,:),1)';
        incoming_conn_I = sum(conns(spikers_I,:),1)';
        G_syn_E(:,t) = G_syn_E(:,t) + parameters.del_G_syn_E*incoming_conn_E; %excitatory conductance
        G_syn_I(:,t) = G_syn_I(:,t) + parameters.del_G_syn_I*incoming_conn_I; %inhibitory conductance
        I_syn(:,t) = G_syn_E(:,t).*(parameters.E_syn_E - V_m(:,t)) + G_syn_I(:,t).*(parameters.E_syn_I - V_m(:,t));
        I_app = I_syn(:,t) + I_theta(:,t);
        %______________________________________
        %Calculate membrane potential using integration method
        V_ss = (I_app + parameters.G_L*parameters.E_L + G_sra(:,t)*parameters.E_K)./(parameters.G_L + G_sra(:,t)); %"steady state" calculation
        taueff = parameters.C_m ./(parameters.G_L + G_sra(:,t)); %timescale for change in the membrane potential
        exp_coeff = (parameters.G_L*(parameters.E_L - V_m(:,t)) + G_sra(:,t).*(parameters.E_K - V_m(:,t)) + I_app)./(parameters.G_L + G_sra(:,t));
        V_m(:,t+1) = V_ss - exp_coeff.*exp(-parameters.dt ./taueff) + randn([parameters.n,1])*(10^(-4)); %the randn portion can be removed if you'd prefer no noise
        V_m(spikers,t+1) = parameters.V_reset; %update those that just spiked to reset
        %______________________________________
        %Update next step conductances
        G_sra(:,t+1) = G_sra(:,t)*exp(-parameters.dt/parameters.tau_sra); %Spike rate adaptation conductance
        %Synaptic conductance updated for each postsynaptic neuron by
        %incoming connection type
        G_syn_E(:,t+1) = G_syn_E(:,t).*exp(-parameters.dt/parameters.tau_syn_E); %excitatory conductance update
        G_syn_I(:,t+1) = G_syn_I(:,t).*exp(-parameters.dt/parameters.tau_syn_I); %inhibitory conductance update
        %______________________________________
        %Update connection strengths
        conns(spikers,:) = parameters.connectivity_gain*conns(spikers,:); %enhance connections of those neurons that just fired
        conns(~spikers,:) = (1/parameters.connectivity_gain)*conns(~spikers,:); %decrease connections of those neurons that did not fire
    end
end