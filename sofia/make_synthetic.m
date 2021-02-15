function [T, U] = make_synthetic(Tsz, R, m, noise)

N = length(Tsz);
U = cell(1, N);

for i=1:N-1
    U{i} = rand(Tsz(i), R);
end
U{N} = zeros(Tsz(N), R);

%TODO: make random pattern
x = 1:1:Tsz(N);
p = linspace(0,1,R+1);
for i=1:R
    amp = sign(rand-0.5)*(1.5*rand + 0.5); % -2 ~ -0.5, 0.5 ~ 2
%     phase_offset = 2*pi*rand; % 0 ~ 2*pi
    phase_offset = 2*pi*p(i);
%     trend = (rand - 0.5) * 0.1; % -0.05 ~ 0.05
    trend = 0;
    level = (rand - 0.5) * 4; % -2 ~ 2
    U{N}(:,i) = amp * sin(2*pi*x/m + phase_offset) + trend*x + level;
end

U = ktensor(U);
[info, ~] = create_problem('Soln', U, 'Noise', noise);
T = info.Data;

end

