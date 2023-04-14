% multivariate t distribution as in Koop book p. 328
function [f]=mstudent(x,m,sigma,v)


k=1; % univariate case only 
c=(pi^(k/2))*gamma(v/2);
cdenom=(v^(v/2))*gamma((v+1)/2);
c=c/cdenom;

f=((det(sigma)).^(-0.5)).*(v+(x-m)'.*inv(sigma).*(x-m))^(-((v+k)/2));
f=f./c;

