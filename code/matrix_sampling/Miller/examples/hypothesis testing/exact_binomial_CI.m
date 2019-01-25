function [lower,upper] = exact_binomial_CI(x,n,confidence,tol)
% Compute an exact confidence interval for p, given an observation X ~ Binomial(n,p).
%
% Input:
%   x = an integer in {0,1,...,n}, representing a single Binomial(n,p) sample.
%   n = a positive integer.
%   confidence = Confidence level to use, e.g. 0.95 for a 95% CI.
%   tol = precision with which to compute the interval endpoints.
% Output:
%   [lower,upper] = lower and upper bounds for the confidence interval.

% Remarks:
% This is a slightly modified version of code by Matt Harrison.

if nargin < 4, tol = 1e-8; end

r = x/n;
a = (1-confidence)/2;

if x == 0
    c1 = 0;
else
    low = 0;
    high = 1;
    c1 = r;
    while high-low > tol
        if 1-binocdf(ceil(x)-1,n,c1) > a
            high = c1;
        else
            low = c1;
        end
        c1 = (low+high)/2;
    end
end

if x == n
    c2 = 1;
else
    low = 0;
    high = 1;
    c2 = r;
    while high-low > tol
        if binocdf(x,n,c2) > a
            low = c2;
        else
            high = c2;
        end
        c2 = (low+high)/2;
    end
end

lower = c1;
upper = c2;