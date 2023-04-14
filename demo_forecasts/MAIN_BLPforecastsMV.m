% ----------------------------------------------------------------------- %
%                                                                         %
%                                                                         %
%          B A Y E S I A N   L O C A L   P R O J E C T I O N              %
%                                                                         %
%  Leonardo N. Ferreira, Silvia Miranda-Agrippino, Giovanni Ricco (2023)  %
%                                                                         %
%                  r e p l i c a t i o n   c o d e                        %
%                                                                         %
%                                                                         %
% Forecast Evaluation                                                     %
% -------------------                                                     %
%                                                                         %                                                                        %
% MODELS:                                                                 %
%              - Local Projections                                        %
%              - Bayesian VAR                                             %
%              - Bayesian Local Projections with RW, VAR and DSGE priors  %
%                                                                         %
% ----------------------------------------------------------------------- %


clear
clc


addpath([pwd '/subroutines/']) %mac
addpath([pwd '/datafolder/'])


%run identifier
runIdentifier='MULTIVARIATE_RW';




%-DATA & MODEL SPECIFICATION----------------------------------------------%



%-LOAD DATA---------------------------------------------------------------%

[tempData,tempText]            =xlsread('DataNew2019.xlsx','FRED_Data_Levels');
dataStructure.data             =tempData(:,2:end);
dataStructure.dates            =x2mdate(tempData(:,1));
dataStructure.varname          =tempText(2,1:end);
dataStructure.description      =tempText(1,1:end);


modelSpec.dataStructure           =dataStructure;


%data characteristics
modelSpec.dataStructure.ragged    =false;


%set estimation sample
modelSpec.beginSample             =datenum(1965,1,1);          


%set forecast sample
modelSpec.beginForecast           =datenum(1990,1,1); 
modelSpec.endForecast             =datenum(2017,12,1);  %last forecast origin


%set forecast type
modelSpec.forecastType            ='recursive'; %recursive/rolling



%initialize BLP prior on presample
modelSpec.presample               =false;%true;

%set BLP pre-sample length (in data points)
modelSpec.presampleLength         =40; %10years


%set BLP prior type
modelSpec.priorType               ='RW'; %BLP only, alternatives: 'VAR', 'AR', 'RW'; 



%set lags & forecast horizon
modelSpec.nVARlags                =5; %number of lags in VAR (& prior)
modelSpec.nBLPlags                =5; %number of lags in BLP
modelSpec.nLPlags                 =5; %number of lags in LP


%forecast horizons
modelSpec.nHorizons               =[1 4 8];






%-HYPERPRIORS SETTINGS----------------------------------------------------%

%random walk prior
hyperPars.isrw                    =true(1,length(modelSpec.dataStructure.varname)); %rw prior

%hyperpriors initial values
%
hyperPars.lambda                  =.4;    %tightness of VAR coeffs prior (the higher lambda the closer to OLS)
hyperPars.lambdaC                 =1e5;   %intercept (you want this to be large)
hyperPars.lambdaP                 =.4;    %tightness of projCoeffs prior (same role as above)
%
hyperPars.miu                     =1;     %sum of coefficients prior (constraint multiplier)
hyperPars.theta                   =2;     %cointegration prior (constraint multiplier) 
hyperPars.alpha                   =2;     %lag decaying coeff for NIW prior


%set hyperpriors options (matches GLP fields); if you want default values
%(when available) set to empty [];
hyperPriorsOptions.hyperpriors    =true;                  %find optimal hyperparameters: NO default option
hyperPriorsOptions.Vc             =1e5;                   %variance of the VAR constant (default=1e6)
hyperPriorsOptions.pos            =find(~hyperPars.isrw); %position of stationary variables
hyperPriorsOptions.MNalpha        =[];                    %lag decaying coeff of NIW prior (default=2)
hyperPriorsOptions.MNpsi          =false;                 %residual variance univariate AR(1) std (default=hyperprior)
hyperPriorsOptions.noc            =false;                 %sum of coefficients prior: NO default option
hyperPriorsOptions.sur            =false;                 %cointegration prior: NO default option
hyperPriorsOptions.Fcast          =false;                 %build forecasts: NO default option
hyperPriorsOptions.hz             =modelSpec.nHorizons;   %max forecast horizon: NO default option
hyperPriorsOptions.mcmc           =false;                 %run metropolis-hasting algorithm: NO default option
hyperPriorsOptions.Ndraws         =1200;                  %default=20k
hyperPriorsOptions.Ndrawsdiscard  =200;                   %default=10k
hyperPriorsOptions.MCMCconst      =1;                     %default=1
hyperPriorsOptions.MCMCfcast      =false;                 %store forecast at each MCMC draw (default=true)
hyperPriorsOptions.MCMCstorecoeff =false;                 %store coefficients at each MCMC draw (default=true)
hyperPriorsOptions.initialValues  =hyperPars;             %see above

hyperPriorsOptions.priorType = modelSpec.priorType;

%sampling from paramenters distribution
GibbsOptions.iterations           =1200;
GibbsOptions.burnin               =200;
GibbsOptions.jump                 =1;

hyperPriorsOptions.GibbsOptions   =GibbsOptions;


modelSpec.hyperPriorsOptions      =hyperPriorsOptions;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%




%-------------------------------------------------------------------------%
%      MULTIVARIATE PSEUDO OUT-OF-SAMPLE FORECASTS FOR LP, VAR, BLP       %
%-------------------------------------------------------------------------%

%BVAR and BLP also produce predictive densities at each horizon

%0) RW
FCST_RW     =forecastRandomWalk(modelSpec);


%1) LP
FCST_LP     =forecastLocalProjMV(modelSpec);


%2) VAR
FCST_VAR    =forecastBayesianVARMV(modelSpec);


%3) BLP
FCST_BLP    =forecastBayesianLocalProjMV(modelSpec);





%-SAVE OUTPUT-------------------------------------------------------------%

modelString =[datestr(modelSpec.beginForecast,'YYYY'),'-',datestr(modelSpec.endForecast,'YYYY')];
outFileName =['FCST_',modelString,'_',runIdentifier,'_',modelSpec.forecastType];


save([pwd '/MATfiles/' outFileName])

%-------------------------------------------------------------------------&


%load([pwd '/MATfiles/FCST_1990-2017_MULTIVARIATE_recursive_paper'])


%%
%-------------------------------------------------------------------------%
%                                 CHARTS                                  % 
%-------------------------------------------------------------------------%
%PLOT FORECASTS

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Courier New')
set(0,'DefaultAxesFontWeight','normal')
set(0,'DefaultAxesFontSize', 8)

% Change default text fonts.
set(0,'DefaultTextFontName', 'Courier New')
set(0,'DefaultTextFontSize', 10)


%
Ocolor  =[255 0 0]/255;
    

%utils
n       =size(modelSpec.dataStructure.data,2);
H       =numel(modelSpec.nHorizons);
labels  =modelSpec.dataStructure.varname;
dates   =modelSpec.dataStructure.dates(...
                modelSpec.dataStructure.dates>modelSpec.beginForecast &...
                modelSpec.dataStructure.dates<=modelSpec.endForecast); 
for h=1:max(modelSpec.nHorizons)
   
    dates=[dates; addtodate(dates(end),1,'month')];
end


startF  =find(modelSpec.dataStructure.dates>=modelSpec.beginForecast,1,'first');
endF    =find(modelSpec.dataStructure.dates<=modelSpec.endForecast,1,'last');

                                   
%plots
figure;         pln=1; 

for h=1:H

    for j=1:n

        subplot(H,n,pln)
        
        hold on
        
        p1=plot(dates,squeeze(FCST_RW.forecasts(:,j,h))./100,':','color',[141 161 177]/255);
        p2=plot(dates,squeeze(FCST_LP.forecasts(:,j,h))./100 ,'--','color',[141 161 177]/255);
        p3=plot(dates,squeeze(FCST_VAR.pointf(:,j,h))./100,'-.','color',[141 161 177]/255);
        p4=plot(dates,squeeze(FCST_BLP.pointf(:,j,h))./100,'-','color',[0 50 100]/255);
        
        p5=plot(dates,modelSpec.dataStructure.data(startF:endF+modelSpec.nHorizons(H)-1,j)./100,'-','color',Ocolor);
        
        set(gca,'FontSize',7,'Xtick',dates(1:36:end),'box','on','layer','top','Ygrid','on'); 
        dateaxis('x',10); datetick('keepticks'); xlim([dates(1) dates(end)]);
        title(labels{j},'Fontweight','Normal','Fontsize',8)
        if j==1 && h==1
            legend([p5 p1 p2 p3 p4],{'data','RW','LP','VAR','BLP'},'FontSize',8,'Location','SouthEast','box','off')
        end
        
        pln=pln+1;
    end
end

set(gcf,'PaperUnits','centimeters','PaperSize',[29 18]) %[x y]
set(gcf,'PaperPosition',[-3.5 0 35 18]) %[left bottom width height]
print(gcf,'-dpdf',[pwd '/CHARTS/PointForecasts_' modelString '.pdf']);            
