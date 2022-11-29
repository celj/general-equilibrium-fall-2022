function f = ll(coef,D,X,M)

alpha = coef(1);
beta = coef(2);
nu = zeros(M, 1);
nu(1:100) = coef(3);

Da = D(:,1);
Db = D(:,2);

p0 = 0.5.*ones(M,2);

options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-6,'FunctionTolerance',1e-6);

probs = fsolve('entry',p0,options,X,alpha,beta,nu);

probsA = probs(:,1);
probsB = probs(:,2);

f = - sum((Da.*Db).*log(probsA.*probsB) ...
    + (Da.*(1 - Db)).*log(probsA.*(1 - probsB)) ...
    + ((1 - Da).*Db).*log((1 - probsA).*probsB) ...
    + ((1 - Da).*(1 - Db)).*log((1 - probsA).*(1 - probsB)));

end