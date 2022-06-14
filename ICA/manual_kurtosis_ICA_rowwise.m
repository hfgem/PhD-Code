function [signals, W, kurt, iter_ns] = manual_kurtosis_ICA_rowwise(X, n_c, iter_n, ...
    alpha, conv_cutoff)

%ABOUT: This function performs a manual ICA implementation using gradient
%descent to maximize kurtosis of components. The function will run for
%either the number of iterations or until change in the cost value of <= 
%cost_cutoff is reached, whichever comes first. This implementation
%independently approaches each signal's weights.
%
%INPUTS:
%   X = dataset to be broken into independent components of size [n_s,n_t],
%       where n_s is the number of signals available and n_t the number of
%       timepoints.
%   n_c = number of independent components to pull out of the data
%   iter_n = number of iterations to run descent
%   alpha = weight of how much mixing values can change at each iteration
%   conv_cutoff = cutoff value of convergence at which to stop iterating
%
%OUTPUTS:
%   signals = matrix of size [n_c,n_t] with all predicted
%       signals based on data samples.
%   W = matrix of size [n_c, n_s] which contains the demixing matrix
%   kurt = average kurtosis of identified signals
%   iter_ns = number of iterations per signal to convergence

    %Grab dimensions
    [n_s,~] = size(X);

    %Initialization
    W = randn(n_c,n_s); %Generate random starting mixing matrix
    W = W ./ sum(W.^2,2); %Have all rows lie on the unit circle

    %Storage matrices
    kurt_vals = zeros(n_c,iter_n); %Store kurtosis values in time for each signal
    cost_vals = zeros(n_c,iter_n); %Store cost values for each signal
    W_last = W;
    cost_last = ones(n_c,1);
    iter_ns = ones(n_c,1); %Store the number of iterations per component

    %1. Calculate starting source signal values
    s_i = W_last*X; %#ok<*MINV>
    %2. Calculate the kurtosis of the starting source signals
    kurt_i = kurtosis(s_i,[],2);
    kurt_vals(:,1) = kurt_i;
    %3. Calculate the cost - recall we want kurtosis to be high for each
    %component and for them to be unique from each other
    cost_i = 1 - abs(tanh(kurt_i));
    cost_vals(:,1) = cost_i;
    %4. Loop through kurtosis ICA implementation
    %Run through each signal independently by estimating only its weights
    for w_i = 1:n_c
        for i = 2:iter_n
            W_test = W_last;
            %5. Test the new cost
            s_i = W_test(w_i,:)*X;
            kurt_i = kurtosis(s_i);
            cost_i = 1 - abs(tanh(kurt_i));
            %6. To maximize kurtosis, we set its derivative = 0 to get the 
            %   following delta_w calculation for updating.
            delta_w = sign(kurt_i)*mean(X.*(s_i).^3,2);
            W_test(w_i,:) = W_test(w_i,:) - alpha*delta_w';
            W_test = W_test./sqrt(sum(W_test.^2,2));
            if w_i > 1
                %To prevent convergence to the same maxima, we must decorelate the 
                %outputs. To do so, we subtract the projections of the previously
                %estimated vectors.
                for j = 1:(w_i-1)
                    W_test(w_i,:) = W_test(w_i,:) - (W_test(w_i,:)*(W_test(j,:)'*W_test(j,:)));
                    W_test(w_i,:) = W_test(w_i,:)./sqrt(W_test(w_i,:)*W_test(w_i,:)');
                end
            end
            delta = abs(sqrt((W_test(w_i,:) - W_last(w_i,:)).^2));
            %7. Store resulting cost, kurtosis, and new matrix
            cost_vals(w_i,i) = cost_i;
            cost_last = cost_i;
            kurt_vals(w_i,i) = kurt_i;
            W_last = W_test;
            if delta <= conv_cutoff
                break
            end
        end
        cost_vals(w_i+1,1) = cost_last;
        iter_ns(w_i,1) = i;
    end
    
    W = W_last;
    signals_pre = W*X; %Predicted signals
    
    %Final separation of signals
    for iter_clean = 1:2
        order = randi(n_c,1,n_c);
        for j = 1:n_c
            row_ind = order(j);
            other_rows = setdiff(1:n_c,row_ind);
            signals_pre(other_rows,:) = signals_pre(other_rows,:) - (signals_pre(other_rows,:)*(signals_pre(row_ind,:)'*signals_pre(row_ind,:))); %Remove other components
            signals_pre = signals_pre./sqrt(sum(signals_pre.^2,2)); %Renormalize
        end
    end
    
    %Store final results
    W = W_last;
    signals = signals_pre; %Predicted signals
    kurt = mean(kurt_vals(:,iter_n));
    
end