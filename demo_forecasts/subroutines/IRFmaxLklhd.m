function res=IRFmaxLklhd(y,modelSpec)

% computes standard ML IRFs 
%
% inputs:
% y=[Txn] matrix of data
% modelSpec=structure with model specification
%
% output: structure (responses to the chosen structural shock only)
% irfs=[(nH+1)xn] matrix of responses to structural shock
% irfs_u & irfs_l=[(nH+1)xn] matrices of error bands
%
% miranda 2014 smirandaagrippino@london.edu

%--------------------------------------------------------------------------

nL=modelSpec.nLags; nH=modelSpec.nHorizons;
[T,n]=size(y); modelSpec.modelSize=n; nB=1000;

%build matrix of relevant lagged y
Ylag=NaN(T-nL,n*nL); %[y_{t-1},...,y_{t-p}]';
for j=1:nL
    Ylag(:,n*(j-1)+1:n*j)=y(nL-j+1:end-j,:);
end
nT=size(Ylag,1); Y=y(nL+1:end,:);

B=[ones(nT,1) Ylag]\Y; VARinnov=Y-[ones(nT,1) Ylag]*B;

MLirfs=SirfBuild(B(:),cov(VARinnov),modelSpec);

%bands
bootInnov=blockBootstrap(VARinnov,nB,1); bootIRFs=NaN(n,nH+1,nB);
for bj=1:nB;
    
    %build bootstrap sample
    ub=squeeze(bootInnov(:,:,bj)); Yb=y; Yb(nL+1:end,:)=NaN;
    for t=1:nT
        
        Yblag=flipud(Yb(t:t+nL-1,:))'; Yblag=[1; Yblag(:)];
        Yb(t+nL,:)=Yblag'*B + ub(t,:);
        
    end
    
    %relevant lags
    Yblag=NaN(T-nL,n*nL); %[y_{t-1},...,y_{t-p}]';
    for j=1:nL
        Yblag(:,n*(j-1)+1:n*j)=Yb(nL-j+1:end-j,:);
    end
    Bb=[ones(nT,1) Yblag]\Yb(nL+1:end,:);
    Sb=cov(Yb(nL+1:end,:)-[ones(nT,1) Yblag]*Bb);
    
    %bootstrapped IRF
    bootIRFs(:,:,bj)=SirfBuild(Bb(:),Sb,modelSpec);
    
    if mod(bj,100)==0
        
        fprintf(1,'bootstrap irfs: sample %i of %i \n',bj,nB);
    end
    
end

bootIRFs=sort(bootIRFs,3);

%bands coverage
bandSize=modelSpec.bandsCoverage;

uBound=bandSize+(100-bandSize)/2; uBound=uBound/100;
lBound=(100-bandSize)/2;          lBound=lBound/100;

%load results
res.irfs   =MLirfs';
res.irfs_l =bootIRFs(:,:,round(uBound*nB))';
res.irfs_u =bootIRFs(:,:,round(lBound*nB))';

%-- child functions -------------------------------------------------------

function irf=SirfBuild(vecB,Sigma,modelSpec)
n=modelSpec.modelSize;
nL=modelSpec.nLags; 
nH=modelSpec.nHorizons;
shockS=modelSpec.shockSize;

Bzero=chol(Sigma); Bzero=bsxfun(@rdivide,Bzero,diag(Bzero));
A=varP2var1(vecB,n,nL);
irf=NaN(n,nH+1); irf(:,1)=eye(n)*Bzero'*shockS;
for h=1:nH
    Ah=A^h; irf(:,h+1)=Ah(1:n,1:n)*Bzero'*shockS;
end

%
function A=varP2var1(vecB,n,nL)
%expects vecB=[n*(1+n*nL)x1]
Ai=reshape(vecB,1+n*nL,n); Ai(1,:)=[]; %remove intercept
A=zeros(n*nL,n*nL); %companion form
A(1:n,:)=Ai'; A(n+1:end,1:n*(nL-1))=eye(n*(nL-1));

%
function bsample=blockBootstrap(x,nsim,blength) % from K. Sheppard
[t,n]=size(x);  x=[x;x(1:blength-1,:)];
s=ceil(t/blength); Bs=floor(rand(s,nsim)*t)+1;

idx=NaN(s*blength,nsim); adder=repmat((0:blength-1)',1,nsim); j=1; 
for i=1:blength:t
    idx(i:(i+blength-1),:)=repmat(Bs(j,:),blength,1)+adder;
    j=j+1;
end
idx=idx(1:t,:);
if n>1
    bsample=NaN(t,n,nsim);
    for i=1:n
        y=x(:,i); bsample(:,i,:)=y(idx);
    end
else
    bsample=x(idx);
end

