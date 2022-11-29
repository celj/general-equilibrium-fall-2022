% Carlos Lezama
% Empirical IO
% Fall 2022

clear;
rng(181121);

%% Data

alpha = 2;
beta = 0.2;
M = 500;
Xa = rand(M, 1);
Xb = 0.1 + (1.4-0.1).*rand(M,1);
X = [Xa Xb];
Ua = randn(M, 1);
Ub = randn(M, 1);
nu = zeros(M, 1);
nu(1:100) = 0.9;

%% Simulation

p0 = 0.5*ones(M,2);

options = optimoptions('fsolve','Display','iter','OptimalityTolerance',1e-6,'FunctionTolerance',1e-6);

probs = fsolve('entry',p0,options,X,alpha,beta,nu);

probsA = probs(:,1);
probsB = probs(:,2);

Da = (norminv(probsA) + Ua > 0);
Db = (norminv(probsB) + Ub > 0);

D = [Da Db];

for i = 1:100:401

    disp(probs(i,:))

end

%% Regret

profitsA = Da.*(Xa*beta - alpha*Db + nu + Ua);
profitsB = Db.*(Xb*beta - alpha*Da + nu + Ub);

entryA = Da(Da(:,:) == 1);
entryB = Db(Db(:,:) == 1);

regretA = profitsA(profitsA(:,:) < 0);
regretB = profitsB(profitsB(:,:) < 0);

disp(size(regretA,1));
disp(size(entryA,1));
disp(size(regretA,1)/size(entryA,1));
disp(size(regretB,1));
disp(size(entryB,1));
disp(size(regretB,1)/size(entryB,1));

%% Estimation

options_ll = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','iter','GradObj','off','HessUpdate','bfgs','TolFun',1e-6,'TolX',1e-6,'MaxIter',1e6,'MaxFunEvals',1e6);

coef0 = ones(3,1);

tic
[coef_est,fval_est,exitflag_est,output,grad_estgh,hessian_estgh] = fminunc('ll',coef0,options_ll,D,X,M);
toc

% SE
se_H_mnl = sqrt(diag(inv(hessian_estgh)));

[coef_est se_H_mnl]
