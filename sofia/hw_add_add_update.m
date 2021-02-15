function [l, b, s] = hw_add_add_update(y, l, b, s, f, m)
% Additive Holt-Winters' algorithm.
%
% INPUTS:
%       y       input vector or matrix.
%       l       level component.
%       b       trend component.
%       s       seasonal component.
%       f       smoothing factors.
%       m       seasonal period.
% OUTPUTS:
%       l       updated level component.
%       b       updated trend component.
%       s       updated seasonal component.
%
% REFERENCE: 
%       D.Lee and K.Shin "Robust Factorization of Real-world Tensor Streams 
%       with Patterns, Missing Values, and Outliers", ICDE 2021.
%
% Created by Dongjin Lee on Feb 15, 2021
% Modified by Dongjin Lee on Feb. 15, 2021

alpha = f(1,:);
beta = f(2,:);
gamma = f(3,:);
alphac = 1 - alpha;
betac = 1 - beta;
gammac = 1 - gamma;
y_alpha = alpha .* y;
y_gamma = gamma .* y;

len = size(l, 1);
n = size(y, 1);

for t=len+1:len+n
    l(t,:) = y_alpha(t-len,:) - alpha .* s(t-m,:) + alphac .* (l(t-1,:) + b(t-1,:));
    b(t,:) = beta .* (l(t,:) - l(t-1,:)) + betac .* b(t-1,:);
    s(t,:) = y_gamma(t-len,:) - gamma .* (l(t-1,:) + b(t-1,:)) + gammac .* s(t-m,:);
end
end
