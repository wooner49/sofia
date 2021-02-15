function results = sofia_wrapper(Y, O, Omega, R, m, cycles, forecast_steps, opts)


%% Initialization
Ysz         = size(Y);
N           = ndims(Y);
ntimes      = Ysz(N);
colons      = repmat({':'}, 1, N-1);

if isfield(opts, 'lambda1'), lambda1 = opts.lambda1;
else, lambda1 = 0.001; end

if isfield(opts, 'lambda2'), lambda2 = opts.lambda2;
else, lambda2 = 0.001; end

if isfield(opts, 'lambda3'), lambda3 = opts.lambda3;
else, lambda3 = 10; end

if isfield(opts, 'mu'), mu = opts.mu;
else, mu = 0.1; end

if isfield(opts, 'phi'), phi = opts.phi;
else, phi = 0.01; end

if isfield(opts, 'maxEpoch'), maxEpoch = opts.maxEpoch;
else, maxEpoch = 500; end

if isfield(opts, 'tol'), tol = opts.tol;
else, tol = 1e-4; end

if isfield(opts, 'Omega_orig'), Omega_orig = tensor(opts.Omega_orig);
else, Omega_orig = tenones(Ysz); end


t_end = Ysz(N)-forecast_steps;
assert(t_end >= m*cycles);

Y_corrupt = Y+O;
Y_corrupt = Y_corrupt(colons{:}, 1:t_end);


%% Method 
[U, W, X_hat, Oout, info_method] = sofia(Y_corrupt, Omega(colons{:}, 1:t_end), ...
    R, m, cycles, lambda1, lambda2, lambda3, mu, phi, maxEpoch, tol, false);


%% Completion results

results.nre = zeros(1, t_end);
results.rmse = zeros(1, t_end);

for t=1:t_end
    results.nre(t) = compute_nre(Omega_orig(colons{:}, t).*Y(colons{:}, t), Omega_orig(colons{:}, t).*X_hat(colons{:}, t));
    results.rmse(t) = compute_rmse(Omega_orig(colons{:}, t).*Y(colons{:}, t), Omega_orig(colons{:}, t).*X_hat(colons{:}, t));
end

results.time = 1:t_end;
results.elapsed_update = info_method.elapsed_dynamic;
results.art_update = results.elapsed_update/ntimes;
results.rae = sum(results.nre(~isinf(results.nre)))/ntimes;
results.U = U;
results.W = W;


%% Forecasting results
if forecast_steps > 0
    results.forecast_time = t_end+1:ntimes;
    L = info_method.L;
    B = info_method.B;
    S = info_method.S;
    U{N} = hw_add_add_forecast(L, B, S, m, forecast_steps);
    Xpred = full(ktensor(U));
    
    results.forecast_nre = zeros(1, forecast_steps);
    results.forecast_rmse = zeros(1, forecast_steps);
    for t=t_end+1:ntimes
        nre = compute_nre(Omega_orig(colons{:}, t).*Y(colons{:}, t), Omega_orig(colons{:}, t).*Xpred(colons{:}, t-t_end));
        rmse = compute_rmse(Omega_orig(colons{:}, t).*Y(colons{:}, t), Omega_orig(colons{:}, t).*Xpred(colons{:}, t-t_end));
        results.forecast_nre(t-t_end) = nre;
        results.forecast_rmse(t-t_end) = rmse;
    end
    results.forecast_rae = sum(results.forecast_nre(~isinf(results.forecast_nre)))/forecast_steps;
end

results.Oout = Oout;
results.internal = info_method;

end