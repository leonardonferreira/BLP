% setpriors.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file sets up the default choices for the priors of the BVAR of 
% Giannone, Lenza and Primiceri (2012)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS:
%
% hyperpriors: 0 = no priors on hyperparameters
%              1 = reference priors on hyperparameters (default)
%              [NOTE: hyperpriors on psi calibrated for data expressed in 
%               4 x logs, such as 4 x log(GDP). Thus if interest rate is in 
%               percentage, divide by 100]     
%
% Vc:       prior variance in the MN prior for the coefficients multiplying
%           the constant term (Default: Vc=10e6)
%           
% pos:      position of the variables that enter the VAR in first
%           differences and for which one might want to set the prior mean 
%           on the coefficient on the first own lag in the MN prior and the
%           prior mean of the sum-of-coefficients prior to 0 (instead of 1)
%           (Default: pos=[])
%
% MNpsi:    0 = diagonal elements of the scale matrix of the IW prior on 
%               the covariance of the residuals NOT treated as 
%               hyperparameters (set to the residual variance of an AR(1))
%           1 = diagonal elements of the scale matrix of the IW prior on 
%               the covariance of the residuals treated as 
%               hyperparameters (default)
%
% MNalpha:  0 = Lag-decaying parameter of the MN prior set to 2 and
%               NOT treated as hyperparameter (default)
%           1 = Lag-decaying parameter of the MN prior treated as 
%               hyperparameter
%
% sur:      0 = single-unit-root prior is OFF
%           1 = single-unit-root prior is ON and its std is treated as an
%               hyperparameter (default)
%
% noc:      0 = no-cointegration (sum-of coefficients) prior is OFF
%           1 = no-cointegration (sum-of coefficients) is ON and its std is 
%               treated as an hyperparameter (default)
%
% fcast:    0 = does not generate forecasts at the posterior mode
%           1 = generates forecasts at the posterior mode (default)
%
% hz:       longest horizon at which the code generates forecasts
%           (default: maxhz=8)
%
% mcmc:     0 = does not run the MCMC (default)
%           1 = runs the MCMC after the maximization
%
% Ndraws:   number of draws in the MCMC (default: Ndraws=20000)
%
% Ndrawsdiscard: number of draws initially discarded to allow convergence 
%                in the in the MCMC (default=Ndraws/2)
%
% MCMCconst: scaling constant for the MCMC (should be calibrated to achieve
%            an acceptance rate of approx 25%) (default: MCMCconst=1)
%
% MCMCfcast:   0 = does not generate forecasts when running the MCMC
%              1 = generates forecasts while running the MCMC
%                  (for each draw of the hyperparameters the code takes a 
%                  draw of the VAR coefficients and shocks, and generates 
%                  forecasts at horizons hz) (default).
%
% MCMCstorecoeff:   0 = does not store the MCMC draws of the VAR 
%                       coefficients and residual covariance matrix
%                   1 = stores the MCMC draws of the VAR coefficients and 
%                       residual covariance matrix (default)
%
% Last modified: 07/01/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Main options

optionsNames=fieldnames(hyperPriorsOptions);

for j=1:length(optionsNames)
    
    eval([optionsNames{j},'=hyperPriorsOptions.',optionsNames{j},';']);
        
end

%--------------------------------------------------------------------------
%switches on optimal prior selection

r.setpriors.hyperpriors=double(hyperpriors);

%--------------------------------------------------------------------------
%variance of prior of constant term  (the larger the wider the variance)
if isempty(Vc)
    Vc=10e6;
end
r.setpriors.Vc=Vc;

%--------------------------------------------------------------------------
%position of variables for which RW prior does not apply

r.setpriors.pos=pos;

%--------------------------------------------------------------------------
%lag decaying rate of the minnesota prior (default = 2), enters prior
%formulation as the exponent to the j-th lag, the higher, the faster the
%decay
if isempty(MNalpha)
    mn.alpha=0;
else
    mn.alpha=MNalpha;
end
r.setpriors.MNalpha=mn.alpha;

%--------------------------------------------------------------------------
%scale of IW for residual variance, if not hyperprior, set to residual 
%variance from univariate AR(1)
if isempty(MNpsi)
    mn.psi=1;
else
    mn.psi=MNpsi;
end
r.setpriors.MNpsi=mn.psi;

%--------------------------------------------------------------------------
%cointegration prior: constraint on optimization which imposes one
%cointegrating vector among ALL variables in the system.

r.setpriors.sur=double(sur);

%--------------------------------------------------------------------------
%sum of coefficients prior: constraint on optimization which imposes FOR
%EACH VARIABLE that the sum of the coefficients of the autoregressive
%polynomial is equal to 1. equivalent to imposing a prior on the SS which
%becomes a function of the unconditional mean of the first p observations.
%irrespective of the value of the constant.

r.setpriors.noc=double(noc);

%--------------------------------------------------------------------------
%builds forecasts at posterior mode

r.setpriors.Fcast=double(Fcast);

%--------------------------------------------------------------------------
%default forecast horizons

r.setpriors.hz=hz;

%--------------------------------------------------------------------------
%MCMC to compute empirical distribution of VAR coefficients: useful for
%IRFs and density forecasts

r.setpriors.mcmc=double(mcmc);

%--------------------------------------------------------------------------
%number of total draws for the MCMC
if isempty(Ndraws)
    M=20000;
else
    M=Ndraws;
end
r.setpriors.Ndraws=M;

%--------------------------------------------------------------------------
%burnin draws
if isempty(Ndrawsdiscard)
    N=round(M/2);
else
    N=Ndrawsdiscard;
end
r.setpriors.Ndrawsdiscard=N;

%--------------------------------------------------------------------------
%scaling constant to calibrate to optimize acceptance rate of MCMC 
if isempty(MCMCconst)
    const=1; 
else
    const=MCMCconst;
end
r.setpriors.MCMCconst=const;

%--------------------------------------------------------------------------
%builds and stores forecasts for each draw
if isempty(MCMCfcast)
    MCMCfcast=1; 
end
r.setpriors.MCMCfcast=double(MCMCfcast);

%--------------------------------------------------------------------------
%stores coefficients at each draw
if isempty(MCMCstorecoeff)
    MCMCstorecoeff=1; 
end
r.setpriors.MCMCstorecoeff=double(MCMCstorecoeff);


%% Other options

% parameters of the hyperpriors, if choosen
if hyperpriors==1;
    
    mode.lambda   =.4;     %hyperpriors modes
    mode.lambdaMA =.1;
    mode.miu      =1;
    mode.theta    =1;
    
    sd.lambda     =.2;     %hyperpriors std
    sd.lambdaMA   =.1;
    sd.miu        =1;
    sd.theta      =1;
    
    scalePSI    =0.02^2; %scale and shape of the IG on psi/(d-n-1)    
    
    %lambda, miu and theta are Gamma distributed with scale (k) and shape
    %(theta) parameters defined by the mode and variance chosen above
    %mode=sd=1 means both the mean and variance are equal to 1
    priorcoef.lambda    =GammaCoef(mode.lambda,sd.lambda,0);  % coefficients of hyperpriors
    priorcoef.lambdaMA  =GammaCoef(mode.lambdaMA,sd.lambdaMA,0);
    priorcoef.miu       =GammaCoef(mode.miu,sd.miu,0);
    priorcoef.theta     =GammaCoef(mode.theta,sd.theta,0);
    
    %this I'm guessing is related to the prior variance of the VAR
    %coefficients but not sure really, I'd expect it to be a function of
    %the residuals variance.. 
    priorcoef.alpha.PSI =scalePSI;
    priorcoef.beta.PSI  =scalePSI;
    
else
    priorcoef=[];
end

%bounds for maximization
MIN.lambda   =0.0001;
MAX.lambda   =5;
%
MIN.lambdaMA =0.0001;
MAX.lambdaMA =5;
%
MIN.theta    =0.0001;
MAX.theta    =50;
%
MIN.miu      =0.0001;
MAX.miu      =50;
%
MIN.alpha    =0.1;
MAX.alpha    =5;

