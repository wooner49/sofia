function TS = tensor2stream(T)

N = ndims(T);
ntimes = size(T, N);

TS = cell(1, ntimes);
colons = repmat({':'}, 1, N);
Tsz = size(T);
Tsz = [Tsz(1:end-1), 1, Tsz(end)];
T = reshape(T, Tsz);

for t = 1:ntimes
    TS{t} = T(colons{:}, t);
end

end

