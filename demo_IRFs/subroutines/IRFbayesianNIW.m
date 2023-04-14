
function res=IRFbayesianNIW(modelSpec,hyperPriorsOptions)

% produces IRFs from BVAR using Gibbs Sampler
% the prior for the VAR coefficients is natural conjugate NIW
%
% inputs:
% modelSpec=structure with model specification
% hyperPars=structure with hyperparameters values
% hyperPriorsOptions=structure with hyperpriors options
%
% output: structure (responses to the chosen structural shock only)
% irfs=[(nH+1)xn] matrix of responses to structural shock
% irfs_u & irfs_l=[(nH+1)xn] matrices of error bands (quantiles)
%
% miranda 2014 smirandaagrippino@london.edu

%--------------------------------------------------------------------------
    
%unload input structures
%model basics
nL =modelSpec.nVARlags; 
nH =modelSpec.nHorizons;

%identification scheme
iScheme=modelSpec.identification;

% %variance of the VAR constant
lambdaC=hyperPriorsOptions.initialValues.lambdaC; %very large number

%find RW variables
isRandomWalk=hyperPriorsOptions.initialValues.isrw;

%Gibbs Sampler
nDraws =hyperPriorsOptions.GibbsOptions.iterations; 
nBurn  =hyperPriorsOptions.GibbsOptions.burnin; 
nJump  =hyperPriorsOptions.GibbsOptions.jump;


%unload data
y      =modelSpec.dataStructure.data;

[T,n]=size(y); modelSpec.modelSize=n;


%build matrix of relevant lagged y
Ylag=NaN(T-nL,n*nL); %[y_{t-1},...,y_{t-p}]';
for j=1:nL
    
    Ylag(:,n*(j-1)+1:n*j)=y(nL-j+1:end-j,:);
    
end

nT=size(Ylag,1); y=y(nL+1:end,:);

y_init=mean(y(1:nL,:)); %average of initial observations


%NIW prior
%Sigma~IW(S_init,a_init)
%vecB|Sigma~N(B_init,V_init) V_init=kron(Sigma,Omega) vecB=[n*(1+n*nL)x1]

%--initialize sampler-----------------------------------------------------%

%set prior residual variance (Sigma) using univariate ar(1) residuals
sigmaj=NaN(n,1); tempYl=reshape(Ylag(:,1:n),nT*n,1);
for j=1:n
    
    sigmaj(j)=std( y(:,j)-[ones(nT,1) tempYl(nT*(j-1)+1:nT*j,:)]*...
        ([ones(nT,1) tempYl(nT*(j-1)+1:nT*j,:)]\y(:,j)) );
    
end

%IW prior for VAR residual variance
a_init=n+2;                 %prior dof (E[Sigma_init]=S_init)

%Gaussian prior for VAR coefficients ~N(B_init,V_init) (equations in columns)
B_init=zeros(n*nL+1,n);
B_init(2:n+1,:)=diag(1.*isRandomWalk); %prior mean VAR coefficients
nB=numel(B_init); 

%projection set
YprojSet=[ones(nT,1) Ylag];


% * * * * * * * * * * * get optimal hyperparameters * * * * * * * * * * * %

parsAtMode=maxMLikelihoodVAR(y,YprojSet,B_init,sigmaj.^2,y_init,hyperPriorsOptions);

lambda=parsAtMode.postmax.lambda; %overall tightness of NIW prior

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %


%Priors' variance
Omega_init=inv(blkdiag(1/lambdaC,kron(diag(1:nL).^2,diag(sigmaj.^2))/lambda^2)); %think of it as inv(Xd'Xd); Xd initial dummy

%-------------------------------------------------------------------------%

%irfs
Omega_end=inv(inv(Omega_init)+YprojSet'*YprojSet); %Kadiyala&Karlsson(1997)

B_end=parsAtMode.postmax.betahat;

v=y-YprojSet*B_end; 

switch iScheme
    case 'CHOL'
        irf=SirfBuild(B_end(:),cov(v),modelSpec);
        
    case 'PSVAR'
        [irf,modelSpec]=psvarIrfBuild(B_end(:),v,modelSpec);
end

%update parameters of the IW
S_end=parsAtMode.postmax.sigmahat*(nT+a_init+n+1); %GLP(2014)

a_end=a_init+nT; %posterior degrees of freedom


%--Gibbs sampler----------------------------------------------------------%

nRetainedDraws=(nDraws-nBurn)/nJump; j=1;
StructuralIRFsCollect=NaN(n,nH+1,nRetainedDraws);

for i=1:nDraws
    
    Sigma_end=iwishrnd(S_end,a_end);    %draw from posterior IW

    %posterior for VAR coefficients ~N(B_end,V_end)    
    %variance
    V_end=kron(Sigma_end,Omega_end);
    
    %mean
    B_end=Omega_end*(Omega_init\B_init+YprojSet'*y);
    
    
    %check VAR stability
    isStableDraw=false;

    while ~isStableDraw
        
        vecB=B_end(:)+chol(V_end)'*randn(nB,1); %draw from posterior N
%         isStableDraw=checkVARstability(vecB,modelSpec);
        isStableDraw=true;
        
    end
    
    v=y-YprojSet*reshape(vecB,numel(vecB)/n,n); 

            
    if i>nBurn && mod(i,nJump)==0 %reduces dependence among draws
        
        switch iScheme
            case 'CHOL'
                StructuralIRFsCollect(:,:,j)=SirfBuild(vecB,Sigma_end,modelSpec); j=j+1;
                
            case 'PSVAR'
                StructuralIRFsCollect(:,:,j)=psvarIrfBuild(vecB,v,modelSpec); j=j+1;
        end
        
    end
    
end

StructuralIRFsCollect=sort(StructuralIRFsCollect,3);


%bands coverage
bandSize=modelSpec.bandsCoverage;

uBound=bandSize+(100-bandSize)/2;   uBound=uBound/100;
lBound=(100-bandSize+1)/2;          lBound=lBound/100;


%load results
res.irfs=irf';
res.irfs_l=StructuralIRFsCollect(:,:,round(lBound*nRetainedDraws))';
res.irfs_u=StructuralIRFsCollect(:,:,round(uBound*nRetainedDraws))';



%-- child functions -------------------------------------------------------
function exitFlag=checkVARstability(vecB,modelSpec)
n=modelSpec.modelSize; nL=modelSpec.nVARlags;

A=varP2var1(vecB,n,nL); exitFlag=true;

if max(abs(eig(A)))>1
    
    exitFlag=false;
    
end


%--------------------------------------------------------------------------
%
function irf=SirfBuild(vecB,Sigma,modelSpec)

n=modelSpec.modelSize;
nL=modelSpec.nVARlags; 
nH=modelSpec.nHorizons;
shockS=modelSpec.shockSize;

Bzero=chol(Sigma); Bzero=bsxfun(@rdivide,Bzero,diag(Bzero));

A=varP2var1(vecB,n,nL);
irf=NaN(n,nH+1); irf(:,1)=eye(n)*Bzero'*shockS; %on impact
for h=1:nH
    
    Ah=A^h; irf(:,h+1)=Ah(1:n,1:n)*Bzero'*shockS;
    
end


%--------------------------------------------------------------------------
%
function A=varP2var1(vecB,n,nL)
%expects vecB=[n*(1+n*nL)x1]

Ai=reshape(vecB,1+n*nL,n); Ai(1,:)=[]; %remove intercept
A=zeros(n*nL,n*nL); %companion form
A(1:n,:)=Ai'; A(n+1:end,1:n*(nL-1))=eye(n*(nL-1));


%--------------------------------------------------------------------------
%
function [irf,modelSpec]=psvarIrfBuild(vecB,VARinnovations,modelSpec)
%irf with external instrument identification

n=modelSpec.modelSize;
nL=modelSpec.nVARlags; 
nH=modelSpec.nHorizons;
shockS=modelSpec.shockSize;

PSVAR=psvar(VARinnovations,modelSpec);
Bzero=PSVAR.Bzero';

% if ~isfield(modelSpec,'Bzero')
% 
%     PSVAR=psvar(VARinnovations,modelSpec);
%     Bzero=PSVAR.Bzero';
%     
%     modelSpec.Bzero =Bzero;
% else
%     Bzero =modelSpec.Bzero;
% end

A=varP2var1(vecB,n,nL);
irf=NaN(n,nH+1); irf(:,1)=eye(n)*Bzero'*shockS; %on impact
for h=1:nH
    
    Ah=A^h; irf(:,h+1)=Ah(1:n,1:n)*Bzero'*shockS;
    
end


%--------------------------------------------------------------------------
%
function res=psvar(VARinnovations,modelSpec)
%identification using external instruments

n  =size(VARinnovations,2);
nL =modelSpec.nVARlags;

instrumentSet=modelSpec.instrument;

%find relevant proxy variable
select=ismember(instrumentSet.labels,{modelSpec.selectedInstrument});
proxy=instrumentSet.data(:,select);


%match sample
residdates =modelSpec.dataStructure.dates(1+nL:end);
proxydates =instrumentSet.dates;

lowerT=max(residdates(1),proxydates(1));
upperT=min(residdates(end),proxydates(end));


innovations =VARinnovations(residdates >= lowerT & residdates <=upperT ,:)';
proxy       =proxy(proxydates >= lowerT & proxydates <=upperT );


keepRows    =~isnan(proxy);
proxy       =proxy(keepRows,:);
innovations =innovations(:,keepRows);

%identification
PSVAR=ProxySVARidentification(innovations,find(modelSpec.shockVar),proxy);

%Bzero
Bzero=zeros(n,n); Bzero(:,modelSpec.shockVar)=PSVAR.B/PSVAR.B(modelSpec.shockVar);



%load structure
res.Bzero       =Bzero;
res.Reliability =PSVAR.L;
res.Fstat       =PSVAR.fstat;

