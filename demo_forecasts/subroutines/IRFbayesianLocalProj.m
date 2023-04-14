
function res=IRFbayesianLocalProj(modelSpec,hyperPriorsOptions)

% computes IRFs using bayesian local projection (SM&GR2015)
%
% Y_{t+h}=a_{h}+B_{1}^{h+1}y_{t-1}+...+B_{p}^{h+1}y_{t-p}+u_{t+h}
% IRF(t,h,d_{i})=E(y_{t+h}|u_t=d_{i})-E(y_{t+h}|u_t=0)=B_{1}^{h}d_{i} h=1:H
% 
% inputs:
% y=[Txn] matrix of data
% modelSpec=structure with model specification
% hyperPriorsOptions=structure with hyperpriors options
%
% output: structure (responses to the chosen structural shock only)
% irfs=[(nH+1)xn] matrix of responses to structural shock
% irfs_u & irfs_l=[(nH+1)xn] matrices of error bands (quantiles)
%
% miranda&ricco 2015 smirandaagrippino@london.edu gricco@london.edu

%--------------------------------------------------------------------------

%unpack input structures
%model basics
nL =modelSpec.nVARlags; %lags in VAR prior
nP =modelSpec.nBLPlags; %lags in BLP

nH=modelSpec.nHorizons; 

shockS=modelSpec.shockSize;

%prior specification
priorType = modelSpec.priorType;

%identification scheme
iScheme=modelSpec.identification;

%variance of the VAR constant
lambdaC=hyperPriorsOptions.initialValues.lambdaC; %very large number

%find RW variables
isRandomWalk=hyperPriorsOptions.initialValues.isrw;

%Gibbs Sampler
GibbsOptions = hyperPriorsOptions.GibbsOptions;

nDraws =GibbsOptions.iterations; 
nBurn  =GibbsOptions.burnin; 
nJump  =GibbsOptions.jump;


%unpack data
y     =modelSpec.dataStructure.data;


%problem size
[T,n]=size(y); modelSpec.modelSize=n;



%-VAR PRIORS---------------------------------------------------------------

%build matrix of relevant lagged y
Ylag=NaN(T-nL,n*nL); %[y_{t-1},...,y_{t-p}]';
for j=1:nL
    
    Ylag(:,n*(j-1)+1:n*j)=y(nL-j+1:end-j,:);
    
end
nT=size(Ylag,1);

Y=y(nL+1:end,:); %dependend Y_{t};

Y_init=mean(y(1:nL,:)); %average of initial observations

%B are the projection coefficients @ different horizons (VAR @ h=0 & h=1)
%Sigma is the covariance of the projection residuals (VAR @ h=0 & h=1)
%Sigma~IW(S,a); vecB|Sigma~N(B,V); V=kron(Sigma,Omega); vecB=[n*(1+n*nLv)x1]


%set prior residual variance (Sigma) using univariate AR(1) residuals
sigmaj=NaN(n,1); 
for j=1:n
    
    sigmaj(j)=std( Y(:,j)-[ones(nT,1) Ylag(:,j)]*...
        ([ones(nT,1) Ylag(:,j)]\Y(:,j)) );
    
end

%IW prior for VAR residual variance
S_init=diag(sigmaj.^2);     %prior scale
a_init=n+2;                 %prior dof (E[Sigma_init]=S_init)

%Gaussian prior for VAR coefficients ~N(B_init,V_init) (equations in columns)
B_init=zeros(n*nL+1,n);
B_init(2:n+1,:)=diag(1.*isRandomWalk); %prior mean (coefficients)

nB=numel(B_init);

%projection set
YprojSet = [ones(nT,1) Ylag]; 


% * * * * * * * * * * * get optimal hyperparameters * * * * * * * * * * * %

parsAtMode=maxMLikelihoodVAR(Y,YprojSet,B_init,sigmaj.^2,Y_init,hyperPriorsOptions);

lambda=parsAtMode.postmax.lambda; %overall tightness of NIW prior

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %

%prior variance (coefficients)
Omega_init=inv(blkdiag(1/lambdaC,kron(diag(1:nL).^2,diag(sigmaj.^2))/lambda^2)); %think of it as inv(Xd'Xd); Xd dummy observations


%--VAR POSTERIOR-----------------------------------------------------------

%posterior variance (coefficients)
Omega_end = inv(inv(Omega_init)+YprojSet'*YprojSet); %Kadiyala&Karlsson(1997)

%posterior mean (coefficients)
B_end = parsAtMode.postmax.betahat;

%update parameters of the IW
%posterior scale
S_end=parsAtMode.postmax.sigmahat*(nT+a_init+n+1);

a_end=a_init+nT; %posterior degrees of freedom

%contemporaneous transmission coefficients
switch iScheme
    
    case 'CHOL'
        Bzero=chol(cov(Y-YprojSet*B_end)); 
        Bzero=bsxfun(@rdivide,Bzero,diag(Bzero));
                
    case 'PSVAR'
        PSVAR=psvar(Y-YprojSet*B_end,modelSpec);
        Bzero=PSVAR.Bzero';
        
        Rho  =PSVAR.Reliability;
        Fstat=PSVAR.Fstat;
        
end


%IRF at posterior mean
irfs=nan(n,nH+1);
irfs(:,1)=shockS'*Bzero;                %on impact
irfs(:,2)=shockS'*Bzero*B_end(2:n+1,:); %h=1


%if VAR/AR-based prior: save coefficients for initialization at future horizons
switch priorType
    
    case 'AR' %use AR(1) coefficients to initialize priors on future horizons
        
        BVARcoeffs=zeros(n*nL+1,n);
        for j=1:n

            BVARcoeffs([1 j+1],j)=[ones(nT,1) Ylag(:,j)]\Y(:,j);

        end    
                

    case 'VAR' %use VAR coefficients to initialize priors on future horizons

        %uses pre-sample
        if modelSpec.presample

            l=modelSpec.nVARlags; %VAR lags in presample initialization

            %estimate VAR on presample
            preY=modelSpec.dataStructure.preSdata;

            preT=size(preY,1);

            %build matrix of relevant lagged y
            preYlag=NaN(preT-l,n*l); %[y_{t-1},...,y_{t-p}]';
            for j=1:l

                preYlag(:,n*(j-1)+1:n*j)=preY(l-j+1:end-j,:);

            end

            preY_init=mean(preY(1:l,:)); %average of initial observations

            preY=preY(l+1:end,:);        %dependend Y_{t};

            priorInit=maxMLikelihoodVAR(preY,[ones(preT-l,1) preYlag],B_init(1:n*l+1,:),sigmaj.^2,preY_init,hyperPriorsOptions);

            %BVAR coefficients on presample
            BVARcoeffs=zeros(n*nL+1,n);
            BVARcoeffs(1:n*l+1,:)=priorInit.postmax.betahat;

        else

            %use all data
            BVARcoeffs=B_end;
            constant(1,:)=B_end(1,:);

        end

end

%----------------------------------------------------------%
%remove deterministic component
tempTrendVAR = handleVARtrend(Y,Ylag,B_end(:),modelSpec);

x = tempTrendVAR.detrended;   %detrended data
%----------------------------------------------------------%

%store sampler's output
nRetainedDraws        =(nDraws-nBurn)/nJump; j=1;
StructuralIRFsCollect =NaN(n,nH+1,nRetainedDraws);
BzeroCollect          =nan(n,n,nRetainedDraws);

RhoCollect            =nan(nRetainedDraws,1);
FstatCollect          =nan(nRetainedDraws,1);

%shrinkage over horizons
optimalLambda=NaN(nH,1); optimalLambda(1)=lambda;


%sampling for VAR coefficients distribution
for i=1:nDraws
    
    Sigma_end=iwishrnd(S_end,a_end);    %draw from posterior IW
    
    %posterior for VAR coefficients ~N(B_end,V_end)    
    %variance
    V_end=kron(Sigma_end,Omega_end);
    
    %mean
    B_end=Omega_end*(Omega_init\B_init+YprojSet'*Y);
    
    %check VAR stability
    isStableDraw=false;

    while ~isStableDraw
        
        vecB=B_end(:)+chol(V_end)'*randn(nB,1); %draw from posterior N
%         isStableDraw=checkVARstability(vecB,modelSpec);        
        isStableDraw=true;
        
    end
    
    %sample Bzero
    switch iScheme
        
        case 'CHOL'
            
            Bzero_=chol(cov( Y-YprojSet*reshape(vecB,numel(vecB)/n,n) )); 
            Bzero_=bsxfun(@rdivide,Bzero_,diag(Bzero_));

        case 'PSVAR'
            
            PSVAR_=psvar(Y-YprojSet*reshape(vecB,numel(vecB)/n,n),modelSpec);
            
            Bzero_=PSVAR_.Bzero';
        
            Rho_  =PSVAR_.Reliability;
            Fstat_=PSVAR_.Fstat;
            
    end
  
        
    %compute and save IRFs
    if i>nBurn && mod(i,nJump)==0 %reduces dependence among draws
        
        %save Bzero for sampling @h>1
        BzeroCollect(:,:,j)=Bzero_;
        
        %store reliability and Fstatistic
        RhoCollect(j)   =Rho_;
        FstatCollect(j) =Fstat_;

        %on impact
        StructuralIRFsCollect(:,1,j)=shockS'*Bzero_; 
        
        %h=1
        vecB=reshape(vecB,numel(vecB)/n,n);
        StructuralIRFsCollect(:,2,j)=shockS'*Bzero_*vecB(2:n+1,:); j=j+1;
                                
    end
    
end


%-BLP AT HORIZON > 1-------------------------------------------------------

%loop over horizons
for h=2:nH
    
    
    if mod(h,10)==0 || h==nH

        system(['say horizon ' num2str(h)]);
    end

    %build relevant projection set (shift obs backward to match horizon)
    XhLag=NaN(T-(nP+h),n*nP);
    for j=h:nP+h-1
    
        XhLag(:,n*(j-h)+1:n*(j-h+1)) = x(nP+h-j+1:end-j,:);

    end  
        
    Xh=x(nP+1+h:end,:); nT=size(Xh,1);
    
    XhprojSet=[ones(nT,1) XhLag]; %XhprojSet=XhLag if no constant
    
    Xh_init=mean(x(1:nP,:)); %mean of initial observations (should now be zero)
    
    
    %use univariate local projection to initialize scale (NW corrected)
    gammaU=NaN(n,1);
    
    for k=1:n
        
        projCoeffs=XhprojSet(:,[1 k+1:n:(n*nP+1)])\Xh(:,k);

        %univariate projection residuals
        u=Xh(:,k)-XhprojSet(:,[1 k+1:n:(n*nP+1)])*projCoeffs;

        u=bsxfun(@minus,u,mean(u)); GammaU=(u'*u)/nT; %HAC error (T*)covariance estimator
        
        nwLags=h-1; %u_{t+h|t} is MA(h-1)
        
        %HAC correction: prior scale
        nwWeights=(nwLags+1-(1:nwLags))./(nwLags+1);
        for j=1:nwLags

            gammaj=(u(j+1:nT,:)'*u(1:nT-j,:))/(nT-j);
            
            GammaU=GammaU+nwWeights(j)*(gammaj+gammaj');

        end

        gammaU(k)=sqrt(GammaU); %scalar
        
    end
    
    
    %prior on proj coeffs: mean
    switch priorType
        
        case 'VAR'
            
            %centered on relevant power of VAR coefficients
            Bh_init=setPriorMean_VAR(BVARcoeffs(2:end,:),h,modelSpec);
        
            
        case 'RW'
            
            %centered on Minnesota-type prior with ones on main diagonal
            Bh_init=setPriorMean_RW(modelSpec,isRandomWalk);
    end

    nBh=numel(Bh_init);
    

    % * * * * * * * * * * get optimal hyperparameters * * * * * * * * * * %

    parsAtMode=maxMLikelihoodBLP(Xh,XhprojSet,Bh_init,gammaU.^2,Xh_init,hyperPriorsOptions,nP,h);

    lambdaP=parsAtMode.postmax.lambda;  %overall tightness of NIW prior
    optimalLambda(h)=lambdaP;

    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %

    
    %prior variance
    OmegaH_init=inv(blkdiag(1/lambdaC,kron(eye(nP),diag(gammaU.^2))/(lambdaP^2))); %think of it as inv(Xd'Xd); Xd dummy observations
        
    %posterior variance
    OmegaH_end=inv(inv(OmegaH_init)+XhprojSet'*XhprojSet); %Kadiyala&Karlsson(1997)

    %posterior mean
    B_end=parsAtMode.postmax.betahat;
    
    constant(h,:)=B_end(1,:);
    
    %IRF at posterior mean
    irfs(:,h+1)=shockS'*Bzero*B_end(2:n+1,:);

    
    %posterior scale
    Sh_end=parsAtMode.postmax.sigmahat*(nT+a_init+n+1);
    
    
    %correct for MA in proj residuals: sandwich covariance matrix for the
    %variance of the projection coefficients
    
    %projection residuals
    Uh=Xh-XhprojSet*B_end;
    Uh=bsxfun(@minus,Uh,mean(Uh)); 

    nwWeights=(nwLags+1-(1:nwLags))./(nwLags+1);
    %HAC correction: posterior scale (mean of posterior IW distribution)
    for l=1:nwLags

        Gammal=(Uh(l+1:nT,:)'*Uh(1:nT-l,:))/(nT-l);
        
        Sh_end=Sh_end+nwWeights(l)*(Gammal+Gammal');

    end
    
    
    %coefficients' distribution
    j=1;
    for i=1:nRetainedDraws
    
        %draw for posterior variance (error)
        SigmaH_end=iwishrnd(Sh_end,a_end);
                
        %posterior variance (coefficients)
        GammaH_end=kron(SigmaH_end,OmegaH_end);
        
        %draw from posterior (coefficients)
        vecB=B_end(:)+chol(GammaH_end)'*randn(nBh,1); %draw from posterior N
        vecB=reshape(vecB,nBh/n,n);
        
        %sample Bzero
        Bzero_=BzeroCollect(:,:,j);
            
        %irf
        StructuralIRFsCollect(:,h+1,j)=shockS'*Bzero_*vecB(2:n+1,:); 

        j=j+1;
                        
    end
end
           
StructuralIRFsCollect =sort(StructuralIRFsCollect,3);
RhoCollect            =sort(RhoCollect);
FstatCollect          =sort(FstatCollect);


%bands coverage
bandSize=modelSpec.bandsCoverage;

uBound=bandSize+(100-bandSize)/2; uBound=uBound/100;
lBound=(100-bandSize)/2;          lBound=lBound/100;


%load results
res.irfs   =irfs';
res.irfs_l =StructuralIRFsCollect(:,:,round(lBound*nRetainedDraws))';
res.irfs_u =StructuralIRFsCollect(:,:,round(uBound*nRetainedDraws))';

res.optimalLambda=optimalLambda;
res.constant     =constant;

res.reliability  =[RhoCollect(round(lBound*nRetainedDraws)) Rho RhoCollect(round(uBound*nRetainedDraws))];
res.Fstatistic   =[FstatCollect(round(lBound*nRetainedDraws)) Fstat FstatCollect(round(uBound*nRetainedDraws))];

%--------------------------------------------------------------------------
%-- child functions -------------------------------------------------------
%--------------------------------------------------------------------------

function exitFlag=checkVARstability(vecB,modelSpec)
n=modelSpec.modelSize; nL=modelSpec.nVARlags;

A=varP2var1(vecB,n,nL); exitFlag=true;

if max(abs(eig(A)))>1
    
    exitFlag=false;
    
end


%--------------------------------------------------------------------------
%
function A=varP2var1(vecB,n,nL)
%expects vecB=[n(n*nL+1)x1]

Ai=reshape(vecB,1+n*nL,n); Ai(1,:)=[]; %remove intercept
A=zeros(n*nL,n*nL); %companion form
A(1:n,:)=Ai'; A(n+1:end,1:n*(nL-1))=eye(n*(nL-1));


%--------------------------------------------------------------------------
%
function res=handleVARtrend(Y,Ylag,vecB,modelSpec)
%removes deterministic component from y; expects B=[n(n*nL+1)x1]

nB=numel(vecB); n=modelSpec.modelSize; nL=modelSpec.nVARlags; nH=modelSpec.nHorizons;

T=size(Y,1);

trend=NaN(T+nH+1,n*nL+1); trend(:,1)=1;
trend(1,2:end)=Ylag(1,:);
for j=1:T+nH
    trend(j+1,2:n+1)=trend(j,:)*reshape(vecB,nB/n,n);
    trend(j+1,n+2:end)=trend(j,2:end-n);
end
trend=trend(:,2:n+1);
x=Y-trend(2:T+1,:);

% figure; plotRows=ceil(n/3); pln=1;
% for j=1:n
%     subplot(plotRows,ceil(n/plotRows),pln)
%     
%     plot(Y(:,j)); hold on
%     plot( trend(2:T,j),'--r'); axis tight
%     pln=pln+1;
% end
% set(gcf,'PaperUnits','centimeters','PaperSize',[18 15]) %[x y]
% set(gcf,'PaperPosition',[-1 0 20 15]) %[left bottom width height]
% print(gcf,'-dpdf','VARtrend.pdf'); saveas(gcf,'VARtrend.fig'); 

res.detrended=[zeros(nL,n);x];
res.trend=trend(2:end,:);       %includes horizons



%--------------------------------------------------------------------------
%
function Bh_init=setPriorMean_RW(modelSpec,isRandomWalk)
%initialize prior on Minnesota-type prior with ones on main diagonal 

%unpack basics
n  =modelSpec.modelSize; 
nP =modelSpec.nBLPlags;


Bh_init=zeros(n*nP+1,n);

%populate diagonal
Bh_init(2:n+1,:)=diag(1.*isRandomWalk); %prior mean (coefficients)


%--------------------------------------------------------------------------
%
function Bh_init=setPriorMean_VAR(BVARcoeffs,h,modelSpec)
%initialize prior on powers of VAR posterior mean BVARcoeffs=[(n*nP)xn] (no
%constant)

%unpack basics
n  =modelSpec.modelSize; 
nL =modelSpec.nVARlags;
nP =modelSpec.nBLPlags;

%companion form
Ai=BVARcoeffs; %do not include constant

A=zeros(n*nL,n*nL); 
A(1:n,:)=Ai'; A(n+1:end,1:n*(nL-1))=eye(n*(nL-1));

%power of VAR coefficients
Bh =A^h; 
Bh =[zeros(1,n); Bh(1:n,:)']; %constant centered around zero


if nL > nP
    Bh_init=Bh(1:n*nP+1,:);
else
    Bh_init=[Bh; zeros(n*(nP-nL),n)];
end


%--------------------------------------------------------------------------
%
function res=psvar(VARinnovations,modelSpec)
%identification using external instruments (Stock&Watson 2012, Mertens&Ravn 2013)

n  =size(VARinnovations,2);
nL =modelSpec.nVARlags;

instrumentSet=modelSpec.instrument;

%find relevant proxy variable
select =ismember(instrumentSet.labels,{modelSpec.selectedInstrument});
proxy  =instrumentSet.data(:,select);


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

