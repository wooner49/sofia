function out = biweight(x)
% k = 1;
% ck = 1.52;

k = 2;
ck = 2.52;
% k = 3;
% ck = 4.12;
out = ck .* ones(size(x));
p = abs(x)<=k;
out(p) = ck * (1 - (1 - (x(p)/k).^2).^3);
end

