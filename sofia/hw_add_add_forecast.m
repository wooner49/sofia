function y_forecast = hw_add_add_forecast(l, b, s, m, h)
% Additive Holt-Winters' algorithm.
%
% INPUTS:
%       l           level component.
%       b           trend component.
%       s           seasonal component.
%       m           seasonal period.
%       h           time steps.
% OUTPUTS:
%       y_forecast  future vector or matrix.
%
% REFERENCE: 
%       D.Lee and K.Shin "Robust Factorization of Real-world Tensor Streams 
%       with Patterns, Missing Values, and Outliers", ICDE 2021.
%
% Created by Dongjin Lee on Feb 15, 2021
% Modified by Dongjin Lee on Feb. 15, 2021

k = size(l, 2);
y_forecast = zeros(h, k);
for t=1:h
    y_forecast(t,:) = l(end,:) + t .* b(end,:) + s(end - m + 1 + mod(t-1, m),:);
end

end

