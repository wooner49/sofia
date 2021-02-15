% compute normalized relative error
function nre = compute_nre(T_true, T_hat)
nre = norm(T_true - T_hat) / norm(T_true);
end

