% Use the exact sampler to estimate the p-value for Galaskiewicz's CEOs and clubs dataset
% using Roberts & Stone (1990)'s $\bar S^2$ statistic.
fprintf('\n===== Galaskiewicz''s CEOs and clubs =====\n');

% number of samples to use
n_samples = 10000;
% confidence (for exact confidence interval)
confidence = 0.95;

matrix_type = 0; % 0:binary, 1:nonnegative integer
A = [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0
     0 0 0 1 0 1 0 0 0 0 1 0 0 1 1 1 1 0 0 1 0 0 1 0 1 1
     1 1 1 1 1 1 1 0 0 1 1 0 1 1 1 1 1 0 1 1 1 1 1 1 1 1
     1 0 0 0 0 0 1 1 1 0 0 1 1 1 0 0 0 1 1 1 1 0 0 1 0 0
     0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0
     0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 1 0 0 0 0
     1 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 0 1 1 0 0 0 0 0 0 0
     0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0
     0 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 1
     0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0 0 1 0
     0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
     0 0 0 1 0 0 0 0 0 0 0 0 0 1 1 0 1 1 0 1 0 1 1 1 0 0];
p = sum(A,2)';
q = sum(A,1);
assert(all(p==[3 11 22 12 3 4 4 4 6 3 4 5 5 3 9]));
assert(all(q==[3 3 2 3 3 3 4 3 4 2 3 2 4 7 5 5 6 5 5 5 3 3 4 5 3 3]));

% Compute the test statistic on the input data
n = length(q);
C = A'*A;
S = sum(sum(triu(C.^2, 1)))/nchoosek(n,2);

% Generate lookup table for binary matrices with row sums p and column sums q
fprintf('Counting...\n');
count(p,q,matrix_type);

% Draw i.i.d. samples from the uniform distribution
fprintf('Sampling %d matrices...\n',n_samples);
X = sample(p,q,n_samples);

% Compute the test statistic on the samples
statistics = zeros(n_samples,1);
for i = 1:n_samples
    C = X(:,:,i)'*X(:,:,i);
    statistics(i) = sum(sum(triu(C.^2, 1)))/nchoosek(n,2);
end

% Compute the p-value
p_value = mean(statistics >= S);
fprintf('Estimated p-value = %f\n',p_value);

% Compute an exact confidence interval for the p-value
[lower,upper] = exact_binomial_CI(n_samples*p_value, n_samples, confidence, 1e-12);
fprintf('Exact %f%% confidence interval: [%f, %f]\n',100*confidence, lower,upper);
fprintf('(%f%% of the time, this interval will contain the true p-value.)\n',100*confidence);











