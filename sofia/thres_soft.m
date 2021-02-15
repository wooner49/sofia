function out = thres_soft(x,v)
x_abs = abs(x);
diff = x_abs - v;
out = sign(x) .* diff .* (diff > 0);
end

