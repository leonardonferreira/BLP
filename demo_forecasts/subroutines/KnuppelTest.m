function [stat,pval] = KnuppelTest(x,lags,prewhite)
% moment test for time series data, using Newey-West covariance matrix

% Inputs:
% x: vector series of PITs (not of S-PITs)
% lags: scalar, lags used for covariance estimation. if ==-1, lag length is chosen
% automatically. Default is lags = -1;
% prewhite: scalar, if == 1, prewhitening used for covariance estimation. Default
% is prewhite = 0.

% Outputs:
% stat: The value of the test statistic KnuppelTest alpha0_1234
% pval: The p-value for KnuppelTest alpha0_1234



s_pit = sqrt(12) * (x - 0.5);
z = [s_pit (s_pit.^2-1)  (s_pit.^3) (s_pit.^4-1.8)];

if nargin < 3;
    prewhite = 0;
    if nargin < 2;
        lags = -1;
    end
end
    
z_orig = z;

z = z_orig(:,[1 3]);
y = sum(z) * (size(z,1)^(-0.5));
[phi,~] = nwX(z-ones(size(z,1),1)*mean(z),prewhite,lags);       
stat_odd = y*inv(phi)*y';

z = z_orig(:,[2 4]);
y = sum(z) * (size(z,1)^(-0.5));
[phi,~] = nwX(z-ones(size(z,1),1)*mean(z),prewhite,lags);       
stat_even = y*inv(phi)*y';

df = 4;
stat = stat_odd + stat_even;
pval = 1 - chi2cdf( stat , df );
        


% Procedure from Bai&Ng, JBES(2005), translated from GAUSS
% u is not demeaned!
function [hac,k] = nwX(u,prewhite,k)

% k is number of lags to be used
% if k is not provided or equals -1, it is determined automatically (Andrews 1991)
% if prewhite == 0 and k == 0, vcv-matrix equals standard vcv-matrix known from OLS 

n = rows(u);
nreg = cols(u);
rho = zeros(nreg,1);
sigma = zeros(nreg,1);
d = zeros(nreg,nreg);
beta = zeros(nreg,nreg);

% do a VAR(1) prewhitening
if prewhite == 1;
    v = zeros(n-1,nreg);
    reg = u(1:n-1,:);
    i = 1;
    while i<=nreg;
        beta(:,i) = lscov( reg , u(2:n,i) );
        v(:,i) = u(2:n,i) - reg * beta(:,i);
        i = i + 1;
    end;
else
    v = u;
end;

if (nargin < 3)||(k==-1);
    i = 1;
    while i<=nreg;
        rho(i) = lscov( v(1:rows(v)-1,i) , v(2:rows(v),i) );
        r = v(2:rows(v),i) - v(1:rows(v)-1,i) * rho(i);
        sigma(i) = r'*r / rows(v);
        i = i + 1;
    end;


    bot = 0; top = 0;
    i = 1;
    while i<=nreg;
        top = top + 4*(rho(i)^2)*sigma(i)^2 / (((1-rho(i))^6)*(1+rho(i))^2);
        bot = bot + sigma(i)^2 / ((1-rho(i))^4);
        i = i+1;
    end;
    alpha = top/bot;
    k = ceil(1.1447*(alpha*n)^(1/3));
end

if k > n/2; 
    k = n/2;
end;

%disp('truncation lag'); 
%disp(k);

vcv = v'*v / (n-1);
i = 1;
while i<=k;
    % x = i/k;
    w = 1 - i/(k+1);
    cov = v(i+1:rows(v),:)' * v(1:rows(v)-i,:) / (n-1);
    vcv = vcv + w*cov;
    cov = v(1:rows(v)-i,:)' * v(i+1:rows(v),:) / (n-1);
    vcv = vcv + w*cov;
    i = i + 1;
end;
d = inv(eye(nreg)-beta');
hac = d*vcv*d';

