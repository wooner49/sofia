function [U,W,X_hat,O,info] = sofia(Y,Omega,R,m,cycles,lambda1,lambda2,lambda3,mu,phi,maxEpoch,tol,needOutlier)
% SOFIA algorithm.
%
% INPUTS:
%       Y            	corrupted tensor 
%       Omega           logical tensor to speficy observable/missing elements.
%       R               rank.
%       m               seasonal period.
%       cycles          number of cycles to be used for SOFIA_INIT.
%       lambda1         hyperparameter.
%       lambda2         hyperparameter.
%       lambda3         hyperparameter.
%       mu              step size.
%       phi             hyperparameter.
%       maxEpoch        max epoch for SOFIA_INIT.
%       tol             tolerance for SOFIA_INIT.
%       needOutlier     flag.
% OUTPUTS:
%       U            	non-temporal factor matrices.
%       W               temporal factor matrix.
%       X_hat           reconstructed tensor.
%       O               outlier tensor.
%       info            information.
%
% REFERENCE: 
%       D.Lee and K.Shin "Robust Factorization of Real-world Tensor Streams 
%       with Patterns, Missing Values, and Outliers", ICDE 2021.
%
% Created by Dongjin Lee on Feb 15, 2021
% Modified by Dongjin Lee on Feb. 15, 2021

if nargin < 11, maxEpoch = 300; end
if nargin < 12, tol = 1e-3; end
if nargin < 13, needOutlier = false; end

fprintf('sofia: start!\n');
%% Parameter setting
Y       = double(Omega .* Y);
Omega   = logical(double((Omega)));
Ysz     = size(Y);
N       = ndims(Y);
ntimes  = Ysz(N);
colons  = repmat({':'}, 1, N-1);


%% Step1: Initialization
ti = m*cycles;
Y_init = Y(colons{:},1:ti);
Omega_init = Omega(colons{:},1:ti);

t_begin_init = tic();
[U_init,X_hat_init,O_init,~] = sofia_init(Y_init,Omega_init,R,m,lambda1,lambda2,lambda3,maxEpoch,tol);

U = cell(N, 1); % Non-temporal factors
W_init = U_init{N}; % Temporal factor
for n=1:N-1 % Normalize
    weights = sqrt(sum(U_init{n}.^2, 1));
    U{n} = U_init{n} ./ weights;
    W_init = W_init .* weights;
end
t_elapsed_init = toc(t_begin_init);

W = zeros(ntimes,R);
W(1:ti,:) = W_init;

X_hat = zeros(size(Y));
X_hat(colons{:},1:ti) = X_hat_init;

O = tensor();
if needOutlier
    O = zeros(size(Y));
    O(colons{:},1:ti)=O_init;
end


%% Step2: HW fitting
t_begin_hw = tic();
[~,L,B,S,F] = hw_add_add_fit(W_init, m);
t_elapsed_hw = toc(t_begin_hw);


%% Initialize error scale tensor
sigma = 0.1.*ones([Ysz(1),Ysz(2)]);

info.scale_time = ti;
info.scale_dynamic_elapsed = 0;

%% Step3: Dynamic Updates
t_begin_dynamic_update = tic();

for t=ti+1:ntimes
    Yt = Y(colons{:},t);
    Omegat = Omega(colons{:},t);
    
    U{N} = hw_add_add_forecast(L,B,S,m,1);
    Ythat = U{1} * diag(U{N}) * U{2}';
    
    Rt= Yt-Ythat;    
    cRt = huber(Rt./sigma).*sigma; % cleaned-residuals
    sigma = sigma_update(sigma, Rt, Omegat, phi);
    
    cRt = Omegat .* cRt;
    
    G = cell(N,1);
    G{1} = cRt * U{2} * diag(U{N});
    G{2} = cRt' * U{1} * diag(U{N});
    G{3} = (khatrirao(U{1},U{2})' * reshape(cRt,[],1))' + lambda1 * (W(t-1,:)-U{N}) + lambda2 * (W(t-m,:)-U{N});
    
    for n=1:N
        G{n} = G{n} * min(1, mu * sqrt(R)/norm(G{n}, 'fro')); % What?
        U{n} = U{n} + mu * G{n};
    end
    
    for n=1:N-1
        weights = sqrt(sum(U{n}.^2, 1));
        U{n} = U{n} ./ weights;
        U{N} = U{N} .* weights;
    end
    
    [L, B, S] = hw_add_add_update(U{N}, L, B, S, F, m);
    W(t,:) = U{N};
    
    if mod(t,500) == 0 % For scalability test.
        info.scale_time = [info.scale_time, t];
        info.scale_dynamic_elapsed = [info.scale_dynamic_elapsed, toc(t_begin_dynamic_update)];
    end
    
    X_hat(colons{:},t) = double(full(ktensor(U)));
    
    if needOutlier
        O(colons{:},t) = Yt-(Yt_hat+cRt);
    end
end
t_elapsed_dynamic = toc(t_begin_dynamic_update);

U = U(1:N-1);

info.elapsed_init = t_elapsed_init;
info.elapsed_hw = t_elapsed_hw;
info.elapsed_dynamic = t_elapsed_dynamic;
info.elapsed_total = t_elapsed_init + t_elapsed_hw + t_elapsed_dynamic;

info.L=L;
info.B=B;
info.S=S;
info.F=F;

fprintf('sofia: end!\n');
end


function new = sigma_update(old, Rt, Omegat, phi)
rho = biweight(Rt./old); 
old_2 = old.^2;

new = phi * rho .* old_2 + (1-phi) * old_2;
new = sqrt(new);
new = Omegat.*new + (1-Omegat).*old;
end
