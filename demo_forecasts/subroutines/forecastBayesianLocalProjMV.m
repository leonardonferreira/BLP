
function res=forecastBayesianLocalProjMV(modelSpec)

% produces BLP-based multivariate OOS forecasts for selected horizons
% BLP as in Miranda-Agrippino & Ricco (2016)
%
% inputs:
% modelSpec     =structure with model specification
%
% output: structure
% BLPforecastAtPostMode =3D container of forecasted values
% BLPforecastErrors     =3D contained of forecast errors
% BLPdensityForecast    =4D container of density forecasts
%
%
% miranda 2019 silvia.miranda-agrippino@northwestern.edu

%--------------------------------------------------------------------------


%unpack model structure
Hstore     =modelSpec.nHorizons;
nH         =numel(Hstore);            %max forecast horizons

fcstType   =modelSpec.forecastType;   %recursive/rolling


sampler    =modelSpec.hyperPriorsOptions.GibbsOptions;


%unpack data
data       =modelSpec.dataStructure.data;
dates      =modelSpec.dataStructure.dates;


%unpack forecast dates
startE     =find(dates>=modelSpec.beginSample,1,'first');
startF     =find(dates>=modelSpec.beginForecast,1,'first');
endF       =find(dates<=modelSpec.endForecast,1,'last');


%unpack presample length (BLP prior)
nO         =modelSpec.presampleLength;
%
n          =size(data,2);
nD         =(sampler.iterations-sampler.burnin)/sampler.jump; %number of draws


%collectors
fcst_collect    =nan(endF-startF+max(Hstore),n,nH);
error_collect   =nan(endF-startF+max(Hstore),n,nH);
density_collect =nan(endF-startF+max(Hstore),n,nH,nD);


%COMPUTE FORECASTS & FORECAST ERRORS
    
for t=startF:endF %for each forecast origin

    %find start date for estimation
    switch fcstType

        case 'recursive'
            beginS =startE;         %sample alway starts from same obs

        case 'rolling'
            beginS =startF-180+1;   %always uses 15 years of data

    end


    y     =data(beginS:t,:);


    % .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
    %pre-sample is the nO observations prior to beginS
    if modelSpec.presample
        py =data(beginS-nO+1:beginS,:); %pre-sample data for BLP initialization

        py(isnan(py),:)=[];

        if numel(py) < nO*n %if less than nO points use the first 5 years of data regardless

            py=data; py(any(isnan(py),2),:)=[]; py=py(1:20,:);
        end

        modelSpec.preSdata =py;
    end
    % .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .


    %produce forecasts
    BLPfcst =forecastBLP(y,modelSpec);


    %STORE
    for h=Hstore %for each horizon

        %forecasts
        fcst_collect(t+h-startF,:,Hstore==h)      = BLPfcst.fcstAtMode(:,h);

        %forecast errors
        error_collect(t+h-startF,:,Hstore==h)     = 1/h*(data(t+h,:) - fcst_collect(t+h-startF,:,Hstore==h));

        %predictive densities
        density_collect(t+h-startF,:,Hstore==h,:) = BLPfcst.fcstDensity(:,h,:);

    end

end
    

%load results
res.pointf      =fcst_collect;
res.errors      =error_collect;
res.densityf    =density_collect;




%--------------------------------------------------------------------------
%-- child functions -------------------------------------------------------
%--------------------------------------------------------------------------

function res =forecastBLP(y,modelSpec)

% produces BLP-based density forecasts
% the prior for the BLP coefficients is NIW centered around either AR or RW


%unpack input structures

%model basics
nL      =modelSpec.nVARlags; %lags in VAR prior
nP      =modelSpec.nBLPlags; %lags in BLP -- if lag selection this becomes the max lag allowed
nH      =max(modelSpec.nHorizons);


%prior specification
priorType = modelSpec.priorType;


%hyperpriors settings
hyperPriorsOptions =modelSpec.hyperPriorsOptions;

%variance of the VAR constant
lambdaC=hyperPriorsOptions.initialValues.lambdaC; %very large number

%find RW variables (all variables have been transformed to to stationarity)
isRandomWalk       =0;


%Gibbs Sampler
nDraws =hyperPriorsOptions.GibbsOptions.iterations; 
nBurn  =hyperPriorsOptions.GibbsOptions.burnin; 
nJump  =hyperPriorsOptions.GibbsOptions.jump;


%-------------------------------------------------------------------------&
%-ESTIMATE BLP-------------------------------------------------------------

%set up data & lags
[T,n] =size(y); 

modelSpec.modelSize=n;


%-VAR PRIORS---------------------------------------------------------------

%build matrix of relevant lagged y
Ylag =nan(T-nL,n*nL); %[y_{t-1},...,y_{t-p}]';
for j=1:nL
    
    Ylag(:,n*(j-1)+1:n*j)=y(nL-j+1:end-j,:);
    
end

nT     =size(Ylag,1);
Y      =y(nL+1:end,:); %dependend Y_{t};

Y_init =mean(y(1:nL,:)); %average of initial observations

%B are the projection coefficients @ different horizons (VAR @ h=0 & h=1)
%Sigma is the covariance of the projection residuals (VAR @ h=0 & h=1)
%Sigma~IW(S,a); vecB|Sigma~N(B,V); V=kron(Sigma,Omega); vecB=[n*(1+n*nLv)x1]


%set prior residual variance (Sigma) using univariate AR(1) residuals
sigmaj=nan(n,1); 
for j=1:n
    
    sigmaj(j)=std( Y(:,j)-[ones(nT,1) Ylag(:,j)]*...
        ([ones(nT,1) Ylag(:,j)]\Y(:,j)) );
    
end

%IW prior for VAR residual variance
S_init          =diag(sigmaj.^2);     %prior scale
a_init          =n+2;                 %prior dof (E[Sigma_init]=S_init)

%Gaussian prior for VAR coefficients ~N(B_init,V_init) (equations in columns)
B_init          =zeros(n*nL+1,n);
B_init(2:n+1,:) =diag(1.*isRandomWalk); %prior mean (coefficients)
nB              =numel(B_init);         %total number of coefficients

%projection set
YprojSet        =[ones(nT,1) Ylag]; 


% * * * * * * * * * * * get optimal hyperparameters * * * * * * * * * * * %

parsAtMode =maxMLikelihoodVAR(Y,YprojSet,B_init,sigmaj.^2,Y_init,hyperPriorsOptions);

lambda     =parsAtMode.postmax.lambda; %overall tightness of NIW prior

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %


%prior variance (coefficients)
Omega_init =inv(blkdiag(1/lambdaC,kron(diag(1:nL).^2,diag(sigmaj.^2))/lambda^2)); %think of it as inv(Xd'Xd); Xd dummy observations




%--VAR POSTERIOR-----------------------------------------------------------
%  @h=1 BLP=VAR

%forecast at mode of the posterior distribution of the parameters

%posterior mean (coefficients)
B_end      =parsAtMode.postmax.betahat;

%posterior variance (coefficients)
Omega_end  =inv(inv(Omega_init)+YprojSet'*YprojSet); %Kadiyala&Karlsson(1997)

%update parameters of the IW
S_end      =parsAtMode.postmax.sigmahat*(nT+a_init+n+1); %posterior scale GLP(2014)

a_end      =a_init+nT;                                   %posterior degrees of freedom




%set up output collectors
Yfcst            =nan(n,nH); %forecast at mode of coefficients' posterior

%sampler output
nD               =(nDraws-nBurn)/nJump; j=1; %number of retained draws
forecastsCollect =nan(n,nH,nD);              %predictive densities


%forecast at horizon 1 (same as VAR)
Ylast      =Y(end:-1:end-nL+1,:)';      Ylast=Ylast(:); %most recent obs

Yfcst(:,1) =[1 Ylast']*B_end; 


%sampling for VAR coefficients distribution
for i=1:nDraws
    
    Sigma_end    =iwishrnd(S_end,a_end);    %draw from posterior IW
    
    %posterior for VAR coefficients ~N(B_end,V_end)    
    %variance
    V_end        =kron(Sigma_end,Omega_end);
    
%     %mean
%     B_end        =Omega_end*(Omega_init\B_init+YprojSet'*Y);
    
    
    %draw coefficients
    vecB         =B_end(:)+chol(V_end)'*randn(nB,1); %draw from posterior N

        
    %store forecasts (distribution)
    if i>nBurn && mod(i,nJump)==0 %reduces dependence among draws
        
        forecastsCollect(:,1,j)=[1 Ylast']*reshape(vecB,n*nL+1,n) + mvnrnd(zeros(n,1),Sigma_end);     
        
        j=j+1;                                   
    end    
end




%--SET BLP PRIORS FOR HORIZONS LARGER THAN 1--


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
            preY=modelSpec.preSdata;

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

        end

end



%--------------------------------------------------------------------------
%-BLP AT HORIZON > 1-------------------------------------------------------

x=y;

%loop over horizons
for h=2:nH
    
    
    if mod(h,10)==0 || h==nH

        system(['say horizon ' num2str(h)]);
    end

    
    %build relevant projection set (shift obs backward to match horizon)
    XhLag =nan(T-(nP+h),n*nP);
    for j=h:nP+h-1
    
        XhLag(:,n*(j-h)+1:n*(j-h+1)) = x(nP+h-j+1:end-j,:);

    end  
        
    Xh        =x(nP+1+h:end,:); 
    nT        =size(Xh,1);
    
    XhprojSet =[ones(nT,1) XhLag]; %XhprojSet=XhLag if no constant
    
    Xh_init   =mean(x(1:nP,:)); %mean of initial observations (should now be zero)
    
    
    %use univariate local projection to initialize scale (NW corrected)
    gammaU=nan(n,1);
    
    for k=1:n
        
        projCoeffs =XhprojSet(:,[1 k+1:n:(n*nP+1)])\Xh(:,k);

        %univariate projection residuals
        u          =Xh(:,k)-XhprojSet(:,[1 k+1:n:(n*nP+1)])*projCoeffs;
        u          =bsxfun(@minus,u,mean(u)); 
        
        GammaU     =(u'*u)/nT; %HAC error (T*)covariance estimator
        
        nwLags     =h+1; 
        
        %HAC correction: prior scale
        nwWeights  =(nwLags+1-(1:nwLags))./(nwLags+1);
        for j=1:nwLags

            gammaj =(u(j+1:nT,:)'*u(1:nT-j,:))/(nT-j);        
            GammaU =GammaU+nwWeights(j)*(gammaj+gammaj');

        end

        gammaU(k)=sqrt(GammaU); %scalar
        
    end
    
    
    %prior on proj coeffs: mean
    switch priorType
        
        case 'VAR'
                        
            %centered on relevant power of VAR coefficients
            Bh_init    =setPriorMean_VAR(BVARcoeffs(2:end,:),h,modelSpec);
        
            
        case 'RW'
            
            %centered on Minnesota-type prior with ones on main diagonal
            Bh_init    =setPriorMean_RW(modelSpec,isRandomWalk);
    end

    nBh =numel(Bh_init);
    

    % * * * * * * * * * * get optimal hyperparameters * * * * * * * * * * %

    parsAtMode  =maxMLikelihoodBLP(Xh,XhprojSet,Bh_init,gammaU.^2,Xh_init,hyperPriorsOptions,nP,h);

    lambdaP     =parsAtMode.postmax.lambda;  %overall tightness of NIW prior

    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %

    
    %posterior mean
    B_end       =parsAtMode.postmax.betahat;
    

    %forecast at posterior mode
    Ylast      =Y(end:-1:end-nP+1,:)';      Ylast=Ylast(:); %most recent obs
    
    Yfcst(:,h)  =[1 Ylast']*B_end;

    
%     %correct for MA in proj residuals: sandwich covariance matrix for the
%     %variance of the projection coefficients
%     
%     %projection residuals
%     Uh      =Xh-XhprojSet*B_end;
%     Uh      =bsxfun(@minus,Uh,mean(Uh)); 
% 
%     nwWeights=(nwLags+1-(1:nwLags))./(nwLags+1);
%     %HAC correction: posterior scale (mean of posterior IW distribution)
%     for l=1:nwLags
% 
%         Gammal =(Uh(l+1:nT,:)'*Uh(1:nT-l,:))/(nT-l);       
%         Sh_end =Sh_end+nwWeights(l)*(Gammal+Gammal');
%     end
    

end
           
forecastsCollect =sort(forecastsCollect,3);


%load results
res.fcstAtMode   =Yfcst;
res.fcstDensity  =forecastsCollect;







%--------------------------------------------------------------------------
%-- child functions -------------------------------------------------------
%--------------------------------------------------------------------------


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




