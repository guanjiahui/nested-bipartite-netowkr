% Use the exact sampler to estimate the p-value for Galton's married couples dataset
% using Diaconis & Efron (1985)'s conditional volume test.

examples_to_run = [1,2];

% number of samples to use
n_samples = 1000;
% confidence (for exact confidence interval)
confidence = 0.95;
matrix_type = 1; % 0:binary, 1:nonnegative integer

for example = examples_to_run
    switch example
        case 1
        % Galton's married couples
        % (Francis Galton, Natural Inheritance, Macmillan, 1889. Page 206.)
        %     s  m  t
        A = [12 20 18   % t = tall
             25 51 28   % m = medium
              9 28 14]; % s = short
        
        case 2
        % Same margins, but different entries
        A = [ 8 14 28
             20 61 23
             18 24  9];

        case 3
        % Double the entries of the previous example
        A = [ 8 14 28
             20 61 23
             18 24  9];
        A = 2*A;
    end
    
    fprintf('\n========== Conditional volume test: Example %d ==========\n',example);

    p = sum(A,2)';
    q = sum(A,1);
            
    % Compute the test statistic on the input data
    total = sum(p);
    M = (p/total)'*(q/total); % proportions under the independence hypothesis
    S = total*sum(sum(((A/total - M).^2)./M)); % chi-square statistic
    
    % Compute Pearson's chi-square test
    m = length(p);
    n = length(q);
    d = (m-1)*(n-1); % degrees of freedom
    chi2_p_value = 1-chi2cdf(S,d);
    fprintf('=== Pearson''s chi-square test ===\n');
    fprintf('Observed chi-square value = %f\n',S);
    fprintf('Degrees of freedom = %d\n',d);
    fprintf('Estimated p-value = %g\n',chi2_p_value);
    fprintf('\n');

    fprintf('=== Conditional volume test ===\n');
    % Generate lookup table for binary matrices with row sums p and column sums q
    fprintf('Counting...\n');
    count(p,q,matrix_type);

    % Draw i.i.d. samples from the uniform distribution
    fprintf('Sampling %d matrices...\n',n_samples);
    X = sample(p,q,n_samples);

    % Compute the test statistic on the samples
    statistics = zeros(n_samples,1);
    for i = 1:n_samples
        statistics(i) = total*sum(sum(((X(:,:,i)/total - M).^2)./M)); % chi-square statistic
    end

    % Compute the p-value
    p_value = mean(statistics <= S);
    fprintf('Estimated p-value = %.12f\n',p_value);

    % Compute an exact confidence interval for the p-value
    [lower,upper] = exact_binomial_CI(n_samples*p_value, n_samples, confidence, 1e-12);
    fprintf('Exact %f%% confidence interval: [%f, %f]\n',100*confidence, lower,upper);
    fprintf('(%f%% of the time, this interval will contain the true p-value.)\n',100*confidence);

end






