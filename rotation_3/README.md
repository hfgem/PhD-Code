<script type="text/x-mathjax-config">
MathJax.Hub.Config({
tex2jax: {
inlineMath: [['$','$'], ['\\(','\\)']],
processEscapes: true},
jax: ["input/TeX","input/MathML","input/AsciiMath","output/CommonHTML"],
extensions: ["tex2jax.js","mml2jax.js","asciimath2jax.js","MathMenu.js","MathZoom.js","AssistiveMML.js", "[Contrib]/a11y/accessibility-menu.js"],
TeX: {
extensions: ["AMSmath.js","AMSsymbols.js","noErrors.js","noUndefined.js"],
equationNumbers: {
autoNumber: "AMS"
}
}
});
</script>

# Rotation 3
 This repository contains code written in Matlab during my third rotation. The goal is to simulate replay spike trains and show that orthagonality is possible within a non-spatial and randomized model of neurons.
 
 The model used in this research is a leaky integrate and fire (LIF) model with spike rate adaptation (SRA) (\cite{miller}). The network is made up of neurons placed in randomly generated clusters, with sparse connectivity, and both excitatory and inhibitory connections. The primary governing equations are as follows:
 
 ![Equations](https://github.com/hfgem/PhD-Code/tree/master/rotation_3/images/equations.png)
 
 Where V_m denotes the membrane potential, V_{reset} the reset membrane potential for the LIF implementation, G_{SRA} denotes the spike rate adaptation conductance, Delta G_{SRA} is a step increase for spike rate adaptation, G_{synE} and G_{synI} are the synaptic conductances specific to excitatory connections and inhibitory connections (respectively), Delta G_{synE/I} are the step increases for the synaptic conductances, $C_m$ is the membrane capacitance, $G_L$ is the leak conductance, $E_L$ is the leak reversal potential, $E_K$ is the potassium reversal potential (associated with the spike rate adaptation current), $I_{syn}$ is the synaptic current, $I_{app}$ is an externally applied current, $\tau_{SRA}$ is the spike rate adaptation decay constant, $\tau_{synE/I}$ are the synaptic decay time constants associated with the excitatory and inhibitory conductances respectively, and $E_{synE/I}$ are the reversal potentials for excitatory and inhibitory currents respectively.
 
 
 
 

 ## lif_clean_code.m
 This file contains code blocks to initialize the model parameters, run the model, and analyze the model outputs. It makes use of the lif_sra_calculator.m function to run the model.
 
 ## lif_sra_calculator.m
 This file contains the function lif_sra_calculator which runs a leaky integrate-and-fire model with spike rate adaptation (sra). The inputs are as follows:
 * n = number of neurons in network
 * clusters = number of clusters in network
 * cluster_mat = number of neurons in a cluster
 * conns = binary nxn matrix of connectivity
 * V_m = n x number of timesteps matrix to store membrane potentials of all neurons over time
 * V_reset = the reset membrane potential
 * V_th = the threshold membrane potential
 * G_sra = n x number of timesteps matrix to store spike rate adaptation conductances of all neurons over time
 * del_G_sra = spike rate adaptation conductance step following spike
 * G_syn_I = n x number of timesteps matrix to store presynaptic inhibitory conductances of all neurons over time
 * G_syn_E = n x number of timesteps matrix to store presynaptic excitatory conductances of all neurons over time
 * del_G_syn_E = excitatory synaptic conductance step following spike
 * del_G_syn_I = inhibitory synaptic conductance step following spike
 * I_syn = n x number of timesteps matrix to store synaptic current values for all neurons over time
 * E_syn_E = nx1 vector of the excitatory synaptic reversal potential
 * E_syn_I = nx1 vector of the inhibitory synaptic reversal potential
 * E_K = potassium reversal potential
 * E_L = leak reversal potential
 * G_L = leak conductance
 * C_m = total membrane capacitance
 * t_start = time step to start the simulation
 * t_end = time step to end the simulation
 * dt = time step size
 * tau_sra = spike rate adaptation time constant
 * tau_syn_E = AMPA/NMDA synaptic decay time constant
 * tau_syn_I = GABA synaptic decay time constant
 * connectivity_gain = amount to increase or decrease connectivity by with each spike (more at the range of 1.002-1.005)
 * type = either 'cluster' or 'neuron' - determines whether or not to initiate sequences with a full cluster of neurons, or a random set of neurons from the network
 * reset = whether or not to reset the neurons to baseline at the simulation start time
 * I_indices = indices of inhibitory neurons
 * E_indices = indices of excitatory neurons
 * I_theta_in = externally applied input (mimicking a theta input)
 * seed = which seed to use for the random number generator for the randomized cluster or neuron selection for initial spiking.