function n = fnorm(T)
v = reshape(T, numel(T), 1);
n = norm(v);
end