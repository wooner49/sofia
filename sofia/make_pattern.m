function P = make_pattern(Tsz, ratio)

if ratio == 0
    P = tensor(true(Tsz));
    return;
end

[info, ~] = create_problem('Size', Tsz, 'M', ratio);
P = tensor(logical(double(info.Pattern)));

end

