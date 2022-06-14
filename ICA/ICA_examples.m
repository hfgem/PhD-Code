%ICA Test Code
%Written by Hannah Germaine
%05-30-2022

%% Create Example Data

%Generate mixture data
dt = 0.01;
t_vec = 0:dt:20;
% s_0 = [sin(40*t_vec); cos(5*t_vec)]; %Sources 1
%s_0 = [2*t_vec; cos(5*t_vec)]; %Sources 2
s_0 = [sawtooth(5*t_vec); sin(-8*t_vec); cos(2*t_vec)]; %Sources 3
s_0 = (s_0 - mean(s_0,2))./std(s_0,[],2); %Source norm
num_samples = 10; %Number of samples of mixture
[n_c,n_t] = size(s_0);
A_0 = rand(num_samples,n_c); %Generated mixture matrix
A_0 = A_0./sqrt(sum(A_0.^2,2));
x = A_0*s_0; %Signal samples

%Data preprocessing - whitening data with SVD
A = x'*x;
epsilon = 0.0001;
[V,D,~] = svd(A);
whMat = sqrt(size(x,1) - 1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
x_white = x*whMat;
de_whMat = pinv(whMat); %To de-whiten perform x_white*de_whMat

%Visualize sources and signal
figure;
for i = 1:n_c+2
    if i <= n_c
        subplot(n_c+2,1,i)
        plot(t_vec,s_0(i,:))
        title(sprintf('True Source %i',i))
    elseif i == n_c + 1
        subplot(n_c+2,1,i)
        plot(t_vec,x)
        title('Mixed Output')
    elseif i == n_c + 2
        subplot(n_c+2,1,i)
        plot(t_vec,x_white)
        title('Whitened Output')
    end
end

%% Perform Manual Kurtosis ICA Implementation
%This will use the kurtosis calculation and gradient descent to move
%towards values for a mixing matrix and source signals that have the
%greatest modulus of kurtosis between them. The cost is greatest the closer
%the modulus of kurtosis is to 0, so we define the cost due to kurtosis as
%   cost = 1/[|kurt(y)| + epsilon]
%where epsilon is a very small value to prevent a division by 0 error.  

%Parameter Definition
iter_n = 1000; %Maximum number of iterations
[n_s,n_t] = size(x_white); %Number of dimensions of data
n_c = 2; %Expected number of independent components
alpha = 0.1; %Weight of how much mixing matrix values can change at each iteration
conv_cutoff = 1*10^(-7); %Iteration cutoff when weight vector change is below this value

%If the number of samples is greater than the number of expected
%components, perform multiple runs of the algorithm on different samples,
%and then average the independent components to find the most likely
%signals
    
[signals, W, kurt, iter_ns] = manual_kurtosis_ICA_rowwise(x_white, n_c, iter_n, ...
alpha, conv_cutoff); %#ok<*ASGLU>

% [signals, W, kurt, iter_run] = manual_kurtosis_ICA(x_white, n_c, iter_n, ...
% alpha, conv_cutoff); %#ok<*ASGLU>


%Plot individual signals
figure;
for i = 1:n_c
    subplot(n_c,1,i)
    plot(t_vec,signals(i,:))
end
sgtitle('All Estimated Signals')

%PCA signals for comparison
[coeff,~,~,~,explained,~] = pca(x_white);
exp_tot = cumsum(explained);
num_sig = find(exp_tot > 95,1);
figure;
for i = 1:num_sig
    subplot(num_sig,1,i)
    plot(t_vec,coeff(:,i)')
end
sgtitle('PCA Components')

%% Perform reconstruction ICA procedure using function
Mdl = rica(x_white,n_c,'IterationLimit',100);
func_s = Mdl.TransformWeights;

figure;
for i = 1:n_c
    subplot(n_c,1,i)
    plot(t_vec,func_s(:,i)')
end
sgtitle('Reconstruction ICA Components')

%% Plot Results

%Plot the resulting signals compared to the original
figure;
ax1 = subplot(4,1,1);
plot(t_vec,s_0)
title('Original Signals')
ax2 = subplot(4,1,2);
plot(t_vec,x_white)
title('Mixed Signal (Whitened)')
ax3 = subplot(4,1,3);
plot(t_vec,signals)
title('Manually Predicted Signals')
ax4 = subplot(4,1,4);
plot(t_vec,func_s')
title('rICA Predicted Signals')

%% Create Example Timeseries Data With Puncta of Activity

%Set up simulation time
t_max = 1; %Maximum simulation time (s)
dt = 0.1*10^(-2); %Timestep (s)
t_vec = 0:dt:t_max; %Time vector

%Set up signal components:
%   1. A base sine wave
base_f = 2; %Base wave frequency (Hz)
base_w = 4; %Base wave amplitude
base_wave = base_w*sin(2*pi*base_f*t_vec) + base_w/2; %Base wave
%   2. A higher frequency sine wave added in random puncta
w1_f = 15; %Wave 1 frequency (Hz)
w1_w = 6; %Wave 1 amplitude
num_puncta = 5; %Number of puncta
punct_len = 0.2; %Length of puncta (s)
w1 = zeros(size(t_vec)); %Wave 1
punct_start = randi(length(t_vec),num_puncta);
for i = 1:num_puncta
    t_start = punct_start(i);
    t_end = punct_start(i)+round(punct_len/dt);
    if t_end > length(t_vec)
        t_end = length(t_vec);
    end
    punct_time = t_vec(t_start:t_end);
    punct_bit = w1_w*sin(2*pi*w1_f*punct_time) + w1_w/2;
    w1(1,t_start:t_end) = punct_bit;
end
%   3. A slow frequency cosine wave
w2_f = 1; %Wave 2 frequency (Hz)
w2_w = 4; %Wave 2 amplitude
w2 = w2_w*cos(2*pi*w2_f*t_vec) + w2_w/2; %Wave 2

%Create multiple combined signals with a given mixture matrix
num_samples = 10; %Number of samples of mixture
num_components = 3;
s_0 = [base_wave;w1;w2];
A_0 = zeros(num_samples,num_components); %Mixture matrix
A_0(:,1) = rand(num_samples,1);
A_0(:,2) = 1 - A_0(:,1);
x = A_0*s_0; %Signal samples

%Whiten data with SVD
A = x'*x;
epsilon = 0.0001;
[V,D,~] = svd(A);
whMat = sqrt(size(x,1) - 1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
x_white = x*whMat;
de_whMat = pinv(whMat); %To de-whiten perform x_white*de_whMat

%Visualize components and combined signal
figure;
ax1 = subplot(5,1,1);
plot(t_vec,base_wave)
title('Base Wave')
xlabel('Time (s)')
ax2 = subplot(5,1,2);
plot(t_vec,w1)
title('Wave 1')
xlabel('Time (s)')
ax3 = subplot(5,1,3);
plot(t_vec,w2)
title('Wave 2')
xlabel('Time (s)')
ax4 = subplot(5,1,4);
plot(t_vec,x)
title('Combined Signal')
xlabel('Time (s)')
ax5 = subplot(5,1,5);
plot(t_vec,x_white);
title('Whitened Signal')
linkaxes([ax1,ax2,ax3])

%% Run Kurtosis ICA on Puncta Data

%Parameter Definition
iter_n = 1000; %Maximum number of iterations
n_c = 3; %Expected number of independent components
[n_s,n_t] = size(x); %Number of dimensions of data
alpha = 0.1; %Weight of how much mixing matrix values can change at each iteration
conv_cutoff = 1*10^(-4); %Iteration cutoff when cost below this value

%Manual Implementation
[signals, W, kurt, iter_ns] = manual_kurtosis_ICA_rowwise(x_white, n_c, iter_n, ...
alpha, conv_cutoff); %#ok<*ASGLU>

% [signals, W, kurt, iter_run] = manual_kurtosis_ICA(x_white, n_c, iter_n, ...
% alpha, conv_cutoff); %#ok<*ASGLU>
    
%Plot individual signals
figure;
for i = 1:n_c
    subplot(n_c,1,i)
    plot(t_vec,signals(i,:))
end
sgtitle('All Estimated Signals')