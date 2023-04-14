
function r=maxMLikelihoodBLP(y,x,b,SS,y0,hyperPriorsOptions,nL,h)
%finds optimal hyperparameters values resulting from maximization of
%log-posterior (marginal likelihood) for application to local projection at 
%further ahead horizons. returns a vector of maximizers and an estimate of 
%the model residual variance. combines together routines by Sims and
%builds on GLP(2014)
%
% inputs:
% y  = [Txn] matrix of observables
% x  = [Tx(1+n*nL)] matrix of lagged observables (constant in first column)
% b  = [(1+n*nL)xn] prior mean for BVAR coefficients
% SS = [nx1] vector of residual univariate AR(1) variances
% y0 = [1xn] vector of mean of initial observations
% hyperPriorsOptions = set of options for the priors specification
%
% output: (structure)
% contains optimal hyperpriors values:
% lambda : NIW prior
% miu    : sum of coefficients prior
% theta  : cointegration prior
% Psi    : residual variance (matrix)
% alpha  : lag decaying coefficient of NIW prior
%
%miranda&ricco (2015) smirandaagrippino@london.edu

%add path to subroutines
addpath([cd '/subroutines'])  %MAC
%addpath([cd '\subroutines']) %PC

%set basic parameters
[T,n] =size(y);
lags  =nL;

%initialize priors - GLP(2014) unmodified
%this isn't exactly true, there's no active role for alpha
setpriors;

%set bounds for residual variance
MIN.psi=SS./100;
MAX.psi=SS.*100;

%-starting values for the maximization-------------------------------------

%----------------------------------------------------%
%OVERWRITES SETTINGS IN SETPRIORS
%Gamma hyperprior for lambda as a function of the horizon
% sd.lambda  =(1+h/5)*.02; %increase variance in h

% L=.35; k=.25;
% sd.lambda=L./(1+exp(-k*(h-5))); %*eig.^(h-10) where eig is max VAR eigenvalue
L=.4; k=.3;
sd.lambda=.1+L./(1+exp(-k*(h-12))); %*eig.^(h-10) where eig is max VAR eigenvalue


mode.lambda=.4;%0.6;         %don't want variance of BLP to be smaller than VAR

%update hyperprior's scale and shape parameters
priorcoef.lambda=GammaCoef(mode.lambda,sd.lambda,0);

priorcoef.priorType = hyperPriorsOptions.priorType;
%----------------------------------------------------%


hyperPars=hyperPriorsOptions.initialValues;

lambda0   =mode.lambda;       %std of NIW prior
theta0    =hyperPars.theta;   %std of cointegration prior (constraint multiplier)
miu0      =hyperPars.miu;     %std of sum of coefficients prior (constraint multiplier)
alpha0    =hyperPars.alpha;   %lag-decaying parameter of the MN prior
psi0      =SS;                %residual variance 


%if hyperpriors is switched on, set starting value for the corresponding hyperparameters
%residual variance
if mn.psi==1;
    inpsi=-log((MAX.psi-psi0)./(psi0-MIN.psi));
elseif mn.psi==0;
    inpsi=[];
end
%lag decaying factor of NIW prior
if mn.alpha==1;
    inalpha=-log((MAX.alpha-alpha0)/(alpha0-MIN.alpha));
elseif mn.alpha==0;
    inalpha=[];
end
%cointegration prior
if sur==1;
    intheta=-log((MAX.theta-theta0)/(theta0-MIN.theta));
elseif sur==0;
    intheta=[];
end
%sum of coefficients prior
if noc==1;
    inmiu=-log((MAX.miu-miu0)/(miu0-MIN.miu));
elseif noc==0;
    inmiu=[];
end

%set initial guess for the hyperparameters vector
x0=[-log((MAX.lambda-lambda0)/(lambda0-MIN.lambda));...
    inpsi;intheta;inmiu;inalpha];
%set initial guess for the corresponding inverse Hessian
H0=10*eye(length(x0));          
%--------------------------------------------------------------------------


%-MAXIMIZE LOG POSTERIOR (MARGINAL LIKELIHOOD)-----------------------------
%returns values of the hyperparametrs which maximize the log posterior
%(conditional on data)
%
[~,hyperParsAtMode,~,~,iterCount,~,~] = csminwel('logMLBLP_formin',...
    x0,H0,...               %initial guesses (parameters and inverse Hessian)
    [],1e-16,1000,...       %gradient,convergence criterion,max number of iterations
    y,x,lags,T,n,...        %dependent,lagged dependent,number of lags,T,n
    b,MIN,MAX,...           %prior mean of BVAR coeffs, optimization bounds
    SS,Vc,...               %prior scale of residual variance, variance of BVAR constant
    pos,mn,...              %position of noRW variables, Minnesota prior's alpha and PSI
    sur,noc,...             %cointegration prior (ON/OFF), sum of coefficients prior (ON/OFF)
    y0,hyperpriors,...      %mean of initial p observations, hyperpriors (ON/OFF)
    priorcoef,h);           %parameters of the hyperpriors distributions
%
%--------------------------------------------------------------------------


%-OPTIMRIZATION OUTPUT: get parameters at mode------------------------------
%VAR coefficients and residual variance
[logML,r.postmax.betahat,r.postmax.sigmahat]=logMLBLP_formin(hyperParsAtMode,y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,hyperpriors,priorcoef,h);

r.lags = lags;                      % # lags
r.postmax.itct=iterCount;           % #iteration before reaching maximum
r.postmax.SSar1=SS;                 % residual variance of AR(1) for each variable
r.postmax.logPost=-logML;           % value of the posterior of the hyperparameters at the peak
r.postmax.lambda=MIN.lambda+(MAX.lambda-MIN.lambda)/(1+exp(-hyperParsAtMode(1)));    % std of MN prior at the peak

r.postmax.theta=MAX.theta;
r.postmax.miu=MAX.miu;

if mn.psi==1;
    
    % diagonal elements of the scale matrix of the IW prior on the residual variance
    r.postmax.psi=MIN.psi+(MAX.psi-MIN.psi)./(1+exp(-hyperParsAtMode(2:n+1)));
    
    if sur==1;
        % std of sur prior at the peak
        r.postmax.theta=MIN.theta+(MAX.theta-MIN.theta)/(1+exp(-hyperParsAtMode(n+2)));
        if noc==1;
            % std of noc prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-hyperParsAtMode(n+3)));
        end
        
    elseif sur==0;
        if noc==1;
            % std of sur prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-hyperParsAtMode(n+2)));
        end
    end
    
elseif mn.psi==0;
    r.postmax.psi=SS;
    
    if sur==1;
        % std of sur prior at the peak
        r.postmax.theta=MIN.theta+(MAX.theta-MIN.theta)/(1+exp(-hyperParsAtMode(2)));
        if noc==1;
            % std of sur prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-hyperParsAtMode(3)));
        end
        
    elseif sur==0;
        if noc==1;
            % std of sur prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-hyperParsAtMode(2)));
        end
    end
end

if mn.alpha==0;
    r.postmax.alpha=2;
elseif mn.alpha==1;
    % Lag-decaying parameter of the MN prior
    r.postmax.alpha=MIN.alpha+(MAX.alpha-MIN.alpha)/(1+exp(-hyperParsAtMode(end)));
end
