function [signals, W, kurt, iter_run] = manual_kurtosis_ICA(X, n_c, iter_n, ...
    alpha, conv_cutoff)

%ABOUT: This function performs a manual ICA implementation using gradient
%descent to maximize kurtosis of components. The function will run for
%either the number of iterations or until change in the cost value of <= 
%cost_cutoff is reached, whichever comes first.
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
%   iter_run = number of iterations run - as a sign of early convergence

    %Grab dimensions
    [n_s,n_t] = size(X);

    %Initialization
    W = rand(n_c,n_s); %Generate random starting mixing matrix
    W = W ./ sum(W.^2,2); %Have all rows lie on the unit circle

    %Storage matrices
    kurt_vals = zeros(n_c,iter_n); %Store kurtosis values in time for each signal
    cost_vals = zeros(n_c,iter_n); %Store cost values for each signal
    W_last = W;
    cost_last = ones(n_c,1);

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
    for i = 2:iter_n
        W_test = W_last;
        %5. Test the new cost
        s_i = W_test*X;
        kurt_i = kurtosis(s_i,[],2);
        cost_i = 1 - abs(tanh(kurt_i));
%             delta = mean(abs(cost_i - cost_last));
        %6. To maximize kurtosis, we set its derivative = 0 to get the 
        %   following delta_w calculation for updating.
        delta_w = sign(kurt_i)'.*(X*(s_i').^3)/n_t;
        W_test = W_test - delta_w';
        W_test = W_test./sqrt(sum(W_test.^2,2));
        %To prevent convergence to the same maxima, we must decorelate the 
        %outputs. To do so, we subtract the projections of the other
        %estimated vectors.
        for j = 1:n_c
            other_rows = setdiff(1:n_c,j);
            W_test(other_rows,:) = W_test(other_rows,:) - (W_test(other_rows,:)*(W_test(j,:)'*W_test(j,:))); %Remove other components
            W_test = W_test./sqrt(sum(W_test.^2,2)); %Renormalize
        end
        %7. Store resulting cost, kurtosis, and new matrix
        cost_vals(:,i) = cost_i;
        cost_last = cost_i;
        kurt_vals(:,i) = kurt_i;
        W_last = W_test;
%             if delta <= conv_cutoff
%                 break
%             end
        if mean(cost_i) <= conv_cutoff
            cost_vals(:,iter_n) = cost_i;
            kurt_vals(:,iter_n) = kurt_i;
            break
        end
    end
    
    W = W_last;
    signals_pre = W*X; %Predicted signals
    
%     %Final separation of signals
%     for s_i = 1:n_c-1
%         s_1 = signals_pre(s_i,:);
%         s_1 = s_1./sqrt(sum(s_1.^2));
%         for s_j = s_i+1:n_c
%             s_2 = signals_pre(s_j,:);
%             s_2 = s_2./sqrt(sum(s_2.^2));
%             s_new = s_1 - (s_1*(s_2'*s_2));
%             s_new = s_new./sqrt(s_new*s_new');
%             signals_pre(s_i,:) = s_new;
%         end
%     end    
    
    %Store final results
    W = W_last;
    signals = signals_pre; %Predicted signals
    kurt = mean(kurt_vals(:,iter_n));
    iter_run = i;
end