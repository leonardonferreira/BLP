
function res=IRFlocalProj(modelSpec)

% computes IRFs using local projection method (Jorda 2005)
%
% Y_{t+h}=a_{h}+B_{1}^{h+1}y_{t-1}+...+B_{p}^{h+1}y_{t-p}+u_{t+h}
% IRF(t,h,d_{i})=E(y_{t+h}|u_t=d_{i})-E(y_{t+h}|u_t=0)=B_{1}^{h}d_{i} h=1:H
% 
% inputs:
% modelSpec=structure with model specification
%
% output: structure (responses to the chosen structural shock only)
% irfs=[(nH+1)xn] matrix of responses to structural shock
% irfs_u & irfs_l=[(nH+1)xn] matrices of error bands 
%
% miranda 2014 smirandaagrippino@london.edu

%--------------------------------------------------------------------------

%unpack model structure
nL =modelSpec.nLPlags; 
nH =modelSpec.nHorizons; 

%identification scheme
iScheme=modelSpec.identification;

%shock size
shockS=modelSpec.shockSize;

%bands size
sLevel=abs(norminv((1-modelSpec.bandsCoverage/100)/2,0,1));


%unpack data
y     =modelSpec.dataStructure.data;

[T,n] =size(y);


%structural responses
LocalProjStructuralIRF=NaN(n,nH+1);
%
lpirfUpperBound=NaN(n,nH+1); 
lpirfLowerBound=NaN(n,nH+1);


%compute IRFs
for h=1:nH

    %build relevant projection set (shift obs backward to match horizon)
    YhLag=NaN(T-(nL+h)+1,n*nL);
    for j=h:nL+h-1

        YhLag(:,n*(j-h)+1:n*(j-h+1)) = y(nL+h-j:end-j,:);

    end  
        
    Yh=y(nL+h:end,:); nT=size(Yh,1);
    
    YhprojSet=[ones(nT,1) YhLag];
    
    iYhprojSet= (YhprojSet'*YhprojSet)\eye(n*nL+1);
        
    %local projection
    projCoeffs=YhprojSet\Yh;
    irfCoeffs=projCoeffs(2:n+1,:); %irf coefficients: equations in columns

    %projection residuals
    u=Yh-YhprojSet*projCoeffs;
    
    %identification
    if h==1
        
        SigmaU=cov(u); 
        %Bzero: structuralShock=VARinnovations*B; B=inv(Bzero);

        switch iScheme
            case 'CHOL'

                Bzero=chol(SigmaU); 
                Bzero=bsxfun(@rdivide,Bzero,diag(Bzero));

            case 'PSVAR'

                Bzero=psvar(u,modelSpec)';
        end
        
        LocalProjStructuralIRF(:,1)=shockS'*Bzero; %on impact
        
    end
    
    %save IRFs
    LocalProjStructuralIRF(:,h+1)=shockS'*Bzero*irfCoeffs;

    
    upperB=LocalProjStructuralIRF(:,h+1);
    lowerB=LocalProjStructuralIRF(:,h+1);
    %robust error bands --Hamilton (1994)
    u=bsxfun(@minus,u,mean(u)); 

    nwLags=h+1;
    nwWeights=(nwLags+1-(1:nwLags))./(nwLags+1);
    
    for jj=1:n
        utemp=u(:,jj);
        err=YhprojSet.*utemp;
        SigmaU=err'*err; 
        
        for j=1:nwLags
            Gammaj=(err(j+1:nT,:)'*err(1:nT-j,:));
            SigmaU=SigmaU+nwWeights(j)*(Gammaj+Gammaj');
        end
        
        sigmaBj = iYhprojSet*SigmaU*iYhprojSet; %HAC error (T*)covariance estimator
        sigmaB = sigmaBj(2:n+1,2:n+1);  %remove only first lags
        
        Omega=(shockS'*Bzero)*sigmaB*(shockS'*Bzero)'; 
    
        Omega=Omega.^0.5*sLevel;

        lpirfUpperBound(jj,h+1)=upperB(jj)+Omega;
        lpirfLowerBound(jj,h+1)=lowerB(jj)-Omega;
    end


end


%load results
res.irfs   =LocalProjStructuralIRF';
res.irfs_u =lpirfUpperBound';
res.irfs_l =lpirfLowerBound';



%-- child functions -------------------------------------------------------
%
function Bzero=psvar(VARinnovations,modelSpec)
%identification using external instruments

n  =size(VARinnovations,2);
nL =modelSpec.nLPlags;

instrumentSet=modelSpec.instrument;

%find relevant proxy variable
select=ismember(instrumentSet.labels,modelSpec.selectedInstrument);
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
res=ProxySVARidentification(innovations,find(modelSpec.shockVar),proxy);

%Bzero
Bzero=zeros(n,n); Bzero(:,modelSpec.shockVar)=res.B/res.B(modelSpec.shockVar);
