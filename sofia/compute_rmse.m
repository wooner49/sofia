% compute root mean squared error
function rmse = compute_rmse(T_true, T_hat)
rmse = sqrt(norm(T_true - T_hat).^2 / prod(size(T_true))); %#ok<PSIZE>
end
