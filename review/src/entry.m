function f = entry(p,X,alpha,beta,nu)

pa = p(:,1);
pb = p(:,2);

Xa = X(:,1);
Xb = X(:,2);

argpa = Xa*beta - alpha*pb + nu;
argpb = Xb*beta - alpha*pa + nu;

f = [pa - normcdf(argpa,0,1); pb - normcdf(argpb,0,1)];

end