function [U,X_hat] = sofia_als(Y,Omega,R,m,lambda1,lambda2,init,maxIters,fitchangetol,printitn)
% SOFIA_ALS algorithm.
%
% INPUTS:
%       Y            	input tensor 
%       Omega           logical tensor to speficy observable/missing elements.
%       R               rank.
%       m               seasonal period.
%       lambda1         hyperparameter.
%       lambda2         hyperparameter.
%       init            initial factor matrices.
%       maxIters        max iters for SOFIA_ALS.
%       fitchangetol    tolerance for SOFIA_ALS.
%       printitn        flag.
% OUTPUTS:
%       U            	factor matrices.
%       X_hat           reconstructed tensor.
%
% REFERENCE: 
%       D.Lee and K.Shin "Robust Factorization of Real-world Tensor Streams 
%       with Patterns, Missing Values, and Outliers", ICDE 2021.
%
% Created by Dongjin Lee on Feb 15, 2021
% Modified by Dongjin Lee on Feb. 15, 2021

%% Parameter setting
Y   = Omega .* Y;
Ysz = size(Y);
N   = ndims(Y);

U_init = init;
for n=1:N-1 % Normalize
    weights = sqrt(sum(U_init{n}.^2, 1));
    U_init{n} = U_init{n} ./ weights;
    U_init{N} = U_init{N} .* weights; 
end


%% Set up for iteration
U = U_init;
normY = fnorm(Y);
fit = 1 - fnorm(Omega .* (Y - double(full(ktensor(U)))))/normY;

if printitn > 0, fprintf('\nsofia_als: start!\n'); end


%% Main Loop
U2 = reshape(U{2}, [1, Ysz(2), 1, R]);
U3 = reshape(U{3}, [1, 1, Ysz(3), R]);

for iter = 1:maxIters
    
    fitold = fit;
          
    % Mode 1 (Non-temporal mode)
    temp1 = U2 .* U3;
    temp2 = tenmat(temp1, N+1).data;
    for i=1:Ysz(1)
        Y_slice = Y(i,:,:);
        Omega_slice = Omega(i,:,:);
        temp3 = reshape(sum(Y_slice .* temp1, [1,2,3]), [1,R]);
        temp4 = temp2(:,Omega_slice);
        temp5 = temp4 * temp4';
        U{1}(i,:) = temp3 * pinv(temp5);
    end
    weights = sqrt(sum(U{1}.^2, 1));
    U{1} = U{1} ./ weights;
    U{N} = U{N} .* weights;
    U1 = reshape(U{1}, [Ysz(1), 1, 1, R]);
    U3 = reshape(U{3}, [1, 1, Ysz(3), R]);
    
    
    % Mode 2 (Non-temporal mode)
    temp1 = U1 .* U3;
    temp2 = tenmat(temp1, N+1).data;
    for i=1:Ysz(2)
        Y_slice = Y(:,i,:);
        Omega_slice = Omega(:,i,:);
        temp3 = reshape(sum(Y_slice .* temp1, [1,2,3]), [1,R]);
        temp4 = temp2(:,Omega_slice);
        temp5 = temp4 * temp4';
        U{2}(i,:) = temp3 * pinv(temp5);
    end
    weights = sqrt(sum(U{2}.^2, 1));
    U{2} = U{2} ./ weights;
    U{N} = U{N} .* weights;
    U2 = reshape(U{2}, [1, Ysz(2), 1, R]);
    
    
    % Mode 3 (Temporal mode)
    temp1 = U1 .* U2;
    temp2 = tenmat(temp1, N+1).data;
    for i=1:Ysz(3)
        Y_slice = Y(:,:,i);
        Omega_slice = Omega(:,:,i);
        temp3 = reshape(sum(Y_slice .* temp1, [1,2,3]), [1,R]);
        temp4 = temp2(:,Omega_slice);
        temp5 = temp4 * temp4';
        
        % Temporal smoothness
        if i <= 1
            temp3 = temp3 + lambda1 * U{3}(i+1,:);
            temp5 = temp5 + lambda1 * eye(R);
        elseif i > 1 && i <= Ysz(3)-1
            temp3 = temp3 + lambda1 * (U{3}(i-1,:) + U{3}(i+1,:));
            temp5 = temp5 + 2 * lambda1 * eye(R);
        else
            temp3 = temp3 + lambda1 * U{3}(i-1,:);
            temp5 = temp5 + lambda1 * eye(R);
        end
        
        % Seasonal smoothness
        if i <= m
            temp3 = temp3 + lambda2 * U{3}(i+m,:);
            temp5 = temp5 + lambda2 * eye(R);
        elseif i > m && i <= Ysz(3)-m
            temp3 = temp3 + lambda2 * (U{3}(i-m,:) + U{3}(i+m,:));
            temp5 = temp5 + 2 * lambda2 * eye(R);
        else
            temp3 = temp3 + lambda2 * U{3}(i-m,:);
            temp5 = temp5 + lambda2 * eye(R);
        end
        U{3}(i,:) = temp3 * pinv(temp5);
    end
    U3 = reshape(U{3}, [1, 1, Ysz(3), R]);
    
    %% Stopping Criteria
    X_hat = double(full(ktensor(U)));
    fit = 1 - fnorm(Omega .* (Y - X_hat))/normY;
    fitchange = abs(fitold - fit);
    
    if (iter > 1) && (fitchange < fitchangetol)
        converge = true;
    else
        converge = false;
    end
    
    if (mod(iter,printitn)==0) || ((printitn>0) && (flag==0))
        fprintf(' sofia_als: iter %2d: fitness = %e, fitness-delta = %7.1e\n', ...
            iter, fit, fitchange);
    end
    
    if converge
        if printitn > 0, fprintf('sofia_als: converge!\n'); end
        break;
    end
end

end