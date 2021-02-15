function [y_hat, l, b, s, f] = hw_add_add_fit(y, m)
% Additive Holt-Winters' algorithm.
%
% INPUTS:
%       y       input vector or matrix.
%       m       seasonal period.
% OUTPUTS:
%       y_hat   reconstructed vector or matrix.
%       l       level component.
%       b       trend component.
%       s       seasonal component.
%       f       smoothing factors.
%
% REFERENCE: 
%       D.Lee and K.Shin "Robust Factorization of Real-world Tensor Streams 
%       with Patterns, Missing Values, and Outliers", ICDE 2021.
%
% Created by Dongjin Lee on Feb 15, 2021
% Modified by Dongjin Lee on Feb. 15, 2021

y_hat = zeros(size(y));
l = zeros(size(y));
b = zeros(size(y));
s = zeros(size(y));
k = size(y, 2);
f = zeros(3, k);

for i=1:k
    [y_hat(:,i), l(:,i), b(:,i), s(:,i), f(:,i)] = hw_add_add_fit_internal(y(:,i), m);
end

end

function [y_hat, l, b, s, f] = hw_add_add_fit_internal(y, m)

assert(isvector(y));
assert(length(y) >= 2*m);

is_y_row = true;
if ~isrow(y)
    is_y_row = false;
    y = y';
end


%% Init Parameters
[l0, b0, s0] = hw_add_add_init_values(y, m);
init_alpha = 0.5 / m;
init_beta = 0.1 * init_alpha;
init_gamma = 0.05 * (1 - init_alpha);

x_init = [init_alpha, init_beta, init_gamma, l0, b0, s0];
% [x_init, max_fval] = hw_add_add_grid_search(x_init, y, m, 20);
lower = [0, 0, 0, -Inf, -Inf, -Inf*ones(size(s0))];
% lower = [0, 0, 0, 0, 0, -Inf*ones(size(s0))];
upper = [1, 1, 1, Inf, Inf, Inf*ones(size(s0))];
max_fval = Inf;
%% Function Handle
funchandle = @(x) hw_add_add_sse_fun(x, y, m, max_fval);


%% 
eps = 1e-4;
bound = 100;
loc = x_init <= lower;
ub = upper;
ub(ub == Inf) = bound;
x_init(loc) = lower(loc) + eps .* (ub(loc) - lower(loc));

loc = x_init >= upper;
lb = lower;
lb(lb == -Inf) = -bound;
x_init(loc) = upper(loc) - eps .*(upper(loc) - lb(loc));

%% Optimization
%TODO: 
options = optimoptions('fmincon', 'HessianApproximation', 'bfgs', ...
                                  'Display', 'off', ...
                                  'MaxFunctionEvaluations', 500000, ...
                                  'MaxIterations', 1000, ...
                                  'SubproblemAlgorithm', 'cg', ...
                                  'OptimalityTolerance', 1e-10, ...
                                  'FiniteDifferenceStepSize', 1e-9, ...
                                  'StepTolerance', 1e-10, ...
                                  'ConstraintTolerance', 1e-10, ...
                                  'FunctionTolerance', 1e-10);

[x,fval,exitflag,output] = fmincon(funchandle, x_init, [], [], [], [], ...
                                   lower, upper, [], options);
fval;
% exitflag
% output
                               
%% Output
f = [x(1), x(2), x(3)];
[y_hat, l, b, s] = hw_add_add_predict(x, y, m);

if ~is_y_row
    y_hat = y_hat';
    l = l';
    b = b';
    s = s';
    f = f';
end

end


function [l0, b0, s0] = hw_add_add_init_values(w, m)
assert(length(w) >= 2*m);

l0 = mean(w(1:m:end));
b0 = mean((w(m+1:2*m) - w(1:m))/m);
s0 = w(1:m) - l0;
end


function [x, minfval] = hw_add_add_grid_search(x, y, m, N)

factors = linspace(0, 1, N);
% factors(1) = 1e-4;
% factors(end) = 1-1e-4;
minfval = Inf;
for a = 1:N
    for b = 1:N
        for c = 1:N
            x(1) = factors(a);
            x(2) = factors(b);
            x(3) = factors(c);
            fval = hw_add_add_sse_fun(x, y, m, Inf);
            if (fval < minfval)
                minfval = fval;
                alpha = factors(a);
                beta = factors(b);
                gamma = factors(c);
            end
        end
    end
end

x(1) = alpha;
x(2) = beta;
x(3) = gamma;

end


function f = hw_add_add_sse_fun(x, y, m, max_fval)

len = length(y);
alpha = x(1);
beta = x(2);
gamma = x(3);
alphac = 1 - alpha;
betac = 1 - beta;
gammac = 1 - gamma;
y_alpha = alpha .* y;
y_gamma = gamma .* y;

if alpha * beta == 0
    f = max_fval;
    return;
end
if beta > alpha || gamma > 1 - alpha
    f = max_fval;
    return;
end

l = zeros(1, len);
b = zeros(1, len);
s = zeros(1, len + m - 1);

l(1) = x(4);
b(1) = x(5);
s(1:m) = x(6:end);

for i=2:len
    l(i) = y_alpha(i-1) - alpha * s(i-1) + alphac * (l(i-1) + b(i-1));
    b(i) = beta * (l(i) - l(i-1)) + betac * b(i-1);
    s(i+m-1) = y_gamma(i-1) - gamma * (l(i-1) + b(i-1)) + gammac * s(i-1);
end

f = norm((l + b + s(1:len)) - y)^2;
end



function [y_hat, l, b, s] = hw_add_add_predict(x, y, m)

len = length(y);
alpha = x(1);
beta = x(2);
gamma = x(3);
alphac = 1 - alpha;
betac = 1 - beta;
gammac = 1 - gamma;
y_alpha = alpha .* y;
y_gamma = gamma .* y;

l = zeros(1, len + 1);
b = zeros(1, len + 1);
s = zeros(1, len + m);

l(1) = x(4);
b(1) = x(5);
s(1:m) = x(6:end);

for i=2:len + 1
    l(i) = y_alpha(i-1) - alpha * s(i-1) + alphac * (l(i-1) + b(i-1));
    b(i) = beta * (l(i) - l(i-1)) + betac * b(i-1);
    s(i+m-1) = y_gamma(i-1) - gamma * (l(i-1) + b(i-1)) + gammac * s(i-1);
end

y_hat = l(1:len) + b(1:len) + s(1:len);
l = l(2:end);
b = b(2:end);
s = s(m+1:end);

end
