function O = make_outlier(Tsz, ratio, magnitude)

O = tensor(rand(Tsz) < ratio);
O = magnitude.*O;
A = (tenrand(Tsz)>0.5)*2-1;
O = A .* O;

end

