
function res=forecastBayesianVAR(modelSpec)

% produces BVAR-based univariate OOS forecasts for selected horizons
%
% inputs:
% modelSpec     =structure with model specification
%
% output: structure
% BVARforecastAtPostMode =3D container of forecasted values
% BVARforecastErrors     =3D contained of forecast errors
% BVARdensityForecast    =4D container of density forecasts
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
leveldata  =modelSpec.dataStructure.leveldata;
datafilter =modelSpec.dataStructure.tcodes;


%unpack forecast dates
startE     =find(dates==modelSpec.beginSample);
startF     =find(dates==modelSpec.beginForecast);
endF       =find(dates==modelSpec.endForecast);


%
n          =size(data,2);


%collectors
fcst_collect    =nan(endF-startF+max(Hstore),n,nH);
error_collect   =nan(endF-startF+max(Hstore),n,nH);
pits_collect    =nan(endF-startF+max(Hstore),n,nH);


%COMPUTE FORECASTS & FORECAST ERRORS

for j=1:n %for each variable
    
    for t=startF:endF %for each forecast origin
            
        %find start date for estimation
        switch fcstType

            case 'recursive'
                beginS =startE;      %sample alway starts from same obs

            case 'rolling'
                beginS =t-startF+1;  %length fixed to first available sample

        end
            
              
        y=data(beginS:t,j);

        if numel(y)>120 && ~any(isnan(y)) %at least 10 years of data & no nans

            

            %produce forecasts
            BVARfcst =forecastBVAR(y,modelSpec);
            
            %revert to original data transform
            BVARfcst.fcstAtMode  =revert(BVARfcst.fcstAtMode,leveldata(1:t,j),datafilter(j));
            BVARfcst.fcstDensity =revert(BVARfcst.fcstDensity,leveldata(1:t,j),datafilter(j));
            
            
            
            
            %STORE
            for h=Hstore %for each horizon
                
                %forecasts
                fcst_collect(t+h-startF,j,Hstore==h)      = BVARfcst.fcstAtMode(:,h);

                %forecast errors
                error_collect(t+h-startF,j,Hstore==h)     = leveldata(t+h,j) - fcst_collect(t+h-startF,j,Hstore==h);
                
                %predictive densities
                pits_collect(t+h-startF,j,Hstore==h)      = ksdensity(BVARfcst.fcstDensity(h,:),leveldata(t+h,j),...
                                                                      'function','cdf','kernel','epanechnikov');
            end
        end
    end
end



%load results
res.pointf      =fcst_collect;
res.errors      =error_collect;
res.pits        =pits_collect;












%--------------------------------------------------------------------------
%-- child functions -------------------------------------------------------
%--------------------------------------------------------------------------


function res =forecastBVAR(y,modelSpec)

% produces BVAR-based density forecasts
% the prior for the VAR coefficients is natural conjugate NIW


%unpack input structures

%model basics
nL      =modelSpec.nVARlags;   %if lag selection this becomes the max lag allowed
nH      =max(modelSpec.nHorizons);


%lag selection
if modelSpec.lagSelection
    
    nL =optimalLag(y,nL,'bic');
    
    %update model specification
    modelSpec.nVARlags =nL;
end


%hyperpriors settings
hyperPriorsOptions =modelSpec.hyperPriorsOptions;


%variance of the VAR constant
lambdaC            =hyperPriorsOptions.initialValues.lambdaC; %very large number

%find RW variables (all variables have been transformed to to stationarity)
isRandomWalk       =0;

%Gibbs Sampler
nDraws             =hyperPriorsOptions.GibbsOptions.iterations; 
nBurn              =hyperPriorsOptions.GibbsOptions.burnin; 
nJump              =hyperPriorsOptions.GibbsOptions.jump;




%-ESTIMATE VAR-------------------------------------------------------------


%set up data & lags
[T,n] =size(y); 


%build matrix of relevant lagged y
ylag =nan(T-nL,n*nL); %[y_{t-1},...,y_{t-p}]';
for j=1:nL
    
    ylag(:,n*(j-1)+1:n*j)=y(nL-j+1:end-j,:);
    
end

nT     =size(ylag,1); 
y      =y(nL+1:end,:);

y_init =mean(y(1:nL,:)); %average of initial observations


%NIW prior
%Sigma~IW(S_init,a_init)
%vecB|Sigma~N(B_init,V_init) V_init=kron(Sigma,Omega) vecB=[n*(1+n*nL)x1]

%--initialize sampler-----------------------------------------------------%

%set prior residual variance (Sigma) using univariate ar(1) residuals
sigmaj=nan(n,1); tempYl=reshape(ylag(:,1:n),nT*n,1);
for j=1:n
    
    sigmaj(j)=std( y(:,j)-[ones(nT,1) tempYl(nT*(j-1)+1:nT*j,:)]*...
        ([ones(nT,1) tempYl(nT*(j-1)+1:nT*j,:)]\y(:,j)) );
    
end

%IW prior for VAR residual variance
a_init =n+2;                 %prior dof (E[Sigma_init]=S_init)

%Gaussian prior for VAR coefficients ~N(B_init,V_init) (equations in columns)
B_init          =zeros(n*nL+1,n);
B_init(2:n+1,:) =diag(1.*isRandomWalk); %prior mean VAR coefficients
nB              =numel(B_init);         %total number of coefficients

%projection set
YprojSet        =[ones(nT,1) ylag];


% * * * * * * * * * * * get optimal hyperparameters * * * * * * * * * * * %

parsAtMode =maxMLikelihoodVAR(y,YprojSet,B_init,sigmaj.^2,y_init,hyperPriorsOptions);

lambda     =parsAtMode.postmax.lambda; %overall tightness of NIW prior

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %


%priors variance
Omega_init =inv(blkdiag(1/lambdaC,kron(diag(1:nL).^2,diag(sigmaj.^2))/lambda^2)); %think of it as inv(Xd'Xd); Xd initial dummy

%-------------------------------------------------------------------------%

%forecast at mode of the posterior distribution of the parameters

B_end      =parsAtMode.postmax.betahat;

Yfcst      =fcstBuild(y,B_end(:),modelSpec);



%update posterior
Omega_end  =inv(inv(Omega_init)+YprojSet'*YprojSet); %Kadiyala&Karlsson(1997)

%update parameters of the IW
S_end      =parsAtMode.postmax.sigmahat*(nT+a_init+n+1); %posterior scale GLP(2014)

a_end      =a_init+nT;                                   %posterior degrees of freedom




%--Gibbs sampler----------------------------------------------------------%

nRetainedDraws   =(nDraws-nBurn)/nJump; j=1;
forecastsCollect =nan(n,nH,nRetainedDraws);


for i=1:nDraws
    
    Sigma_end    =iwishrnd(S_end,a_end);    %draw from posterior IW

    %posterior for VAR coefficients ~N(B_end,V_end)    
    %variance
    V_end        =kron(Sigma_end,Omega_end);
    
    %mean
    B_end        =Omega_end*(Omega_init\B_init+YprojSet'*y);
    
    %draw coefficients
    vecB         =B_end(:)+chol(V_end)'*randn(nB,1); %draw from posterior N

            
    if i>nBurn && mod(i,nJump)==0 %reduces dependence among draws
        
        
        forecastsCollect(:,:,j)=fcstBuildSampling(y,vecB,Sigma_end,modelSpec);     j=j+1;                
    end    
end

forecastsCollect =sort(forecastsCollect,3);



%load results
res.fcstAtMode  =Yfcst;
res.fcstDensity =forecastsCollect;




%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


function fcst =fcstBuild(y,vecB,modelSpec)

nL =modelSpec.nVARlags; 
nH =max(modelSpec.nHorizons);


[T,n]=size(y);


%last observation
yT =[1 y(end:-1:end-nL+1,:)'];


%store forecasts
fcst =nan(n,T+nH);
fcst(:,1:T) =y'; 

for h=1:nH
        
    fcst(:,T+h) =yT*vecB;    
    
    yT =[1 fcst(:,T+h:-1:T+h-nL+1)];    
end

%load results
fcst =fcst(:,T+1:end);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function fcst =fcstBuildSampling(y,vecB,Sigma,modelSpec)

nL =modelSpec.nVARlags; 
nH =max(modelSpec.nHorizons);


[T,n] =size(y);


%store forecasts
fcst =nan(n,T+nH);
fcst(:,1:T) =y'; 

for h=1:nH
        
    yT =fcst(:,T+h-1:-1:T+h-nL);    yT=yT(:);
    
    fcst(:,T+h) =[1 yT']*reshape(vecB,n*nL+1,n) + mvnrnd(zeros(n,1),Sigma);    
    
end

%load results
fcst =fcst(:,T+1:end);



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function res =optimalLag(y,nL,criterion)

Ocrit =nan(nL,1);

[T,n] =size(y); 

l     =1;


while l <= nL 
        
    ylag =nan(T-l,n*l); %[y_{t-1},...,y_{t-p}]';
    for j=1:l

        ylag(:,n*(j-1)+1:n*j)=y(l-j+1:end-j,:);

    end

    nT     =size(ylag,1); 
    Y      =y(l+1:end,:);

    YprojSet =[ones(nT,1) ylag];


    switch criterion

        case 'aic'
            Ocrit(l) =aic(Y,YprojSet);

        case 'bic'
            Ocrit(l) =bic(Y,YprojSet);
    end
    
    if l>1 && find(Ocrit==min(Ocrit)) ~= sum(~isnan(Ocrit))
        
        break
    end
    l=l+1;
end

res =find(Ocrit==min(Ocrit));



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function res =revert(modelfcst,data,tcode)

%data at forecast origin: in levels
y     =data(end);
dy    =data(end)-data(end-1);

ly    =log(data(end));
dly   =log(data(end))-log(data(end-1));


%model forecast 
f_old =modelfcst;  %1 x nH x nD

%max forecast horizon
nH    =size(f_old,2);

%number of draws
nD    =size(f_old,3);


%revert transformation
%random walk is treated as a direct forecast
switch tcode
   
    case 1 % level (i.e. no transformation): x(t)
        
        f_new =f_old;
        
        
    case 2 % first difference: x(t)-x(t-1)
        
        f_new        =nan(1,nH+1,nD);
        f_new(:,1,:) =y;
        
        for h=1:nH
           
            f_new(:,h+1,:) =f_new(:,h,:) + f_old(:,h,:); 
        end        
        f_new   =f_new(:,2:end,:);
        
        
    case 3 % second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
        
        f_new        =nan(1,nH+1,nD);
        f_new(:,1,:) =y + dy;

        for h=1:nH
           
            f_new(:,h+1,:) =f_new(:,h,:) + dy*ones(1,nD) + f_old(:,h,:); 
        end        
        f_new   =f_new(:,2:end,:);
        
        
    case 4 % natural log: ln(x)
        
        f_new =exp(f_old);
        
        
    case 5 % first difference of natural log: ln(x)-ln(x-1)
        
        f_new        =nan(1,nH+1,nD);        
        f_new(:,1,:) =ly;
        
        for h=1:nH
           
            f_new(:,h+1,:) =f_new(:,h,:) + f_old(:,h,:); 
        end     
        f_new =exp(f_new(:,2:end,:));
        
        
    case 6 % second difference of natural log: (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
        
        f_new        =nan(1,nH+1,nD);        
        f_new(:,1,:) =ly + dly;
        
        for h=1:nH
           
            f_new(:,h+1,:) =f_new(:,h,:) + dly*ones(1,1,nD) + f_old(:,h,:);
        end
        f_new =exp(f_new(:,2:end,:));
                
end


%load results
res = squeeze(f_new);




