
% ----------------------------------------------------------------------- %
%                                                                         %
%                                                                         %
%          B A Y E S I A N   L O C A L   P R O J E C T I O N              %
%                                                                         %
%  Leonardo N. Ferreira, Silvia Miranda-Agrippino, Giovanni Ricco (2023)  %
%                                                                         %
%                                                                         %
%                  r e p l i c a t i o n   c o d e                        %
% ----------------------------------------------------------------------- %

clear
clc

addpath([pwd '/subroutines/']) %mac


%run identifier
runIdentifier='SUBSAMPLERW';


%-LOAD DATA---------------------------------------------------------------%

[tempData,tempText]      =xlsread('DataNew2019.xlsx','FRED_Data_Levels');
VARdata.data             =tempData(:,2:end);
VARdata.dates            =x2mdate(tempData(:,1));
VARdata.labels           =tempText(2,1:end);
VARdata.description      =tempText(1,1:end);



%-MODEL SPECIFICATION-----------------------------------------------------%

%choose VAR variables
varList                        =VARdata.labels;

%initialize BLP prior on presample
modelSpec.presample            =true;

%set BLP prior type
modelSpec.priorType            ='RW'; %BLP only, alternatives: 'VAR', 'RW'; 

modelSpec.nSubsamples          =25; %loops over rolling 30-years subsamples
%first run: 1965-1995
%last  run: 1989-2019


%set VAR specification, lags, IRFs horizon
modelSpec.nVARlags             =5; %number of lags in VAR (& prior)
modelSpec.nBLPlags             =5; %number of lags in BLP
modelSpec.nLPlags              =5; %number of lags in LP

modelSpec.nHorizons            =20;
modelSpec.bandsCoverage        =90;


%choose identification scheme
modelSpec.identification       ='CHOL'; %alternatives: 'PSVAR', 'CHOL';


%declare shock variable & shock size
shockVar                       ='FFR';       %1-year rate
shockSize                      =1;           %percentage points


%-HYPERPRIORS SETTINGS----------------------------------------------------%

%random walk prior
hyperPars.isrw                   =true(1,numel(varList)); %rw prior

%hyperpriors initial values
%
hyperPars.lambda                 =.4;    %tightness of VAR coeffs prior (the higher lambda the closer to OLS)
hyperPars.lambdaC                =1e5;   %intercept (you want this to be large)
hyperPars.lambdaP                =.4;    %tightness of projCoeffs prior (same role as above)
%
hyperPars.miu                    =1;     %sum of coefficients prior (constraint multiplier)
hyperPars.theta                  =2;     %cointegration prior (constraint multiplier) 
hyperPars.alpha                  =2;     %lag decaying coeff for NIW prior


%set hyperpriors options (matches GLP fields); if you want default values
%(when available) set to empty [];
hyperPriorsOptions.hyperpriors   =true;                  %find optimal hyperparameters: NO default option
hyperPriorsOptions.Vc            =1e5;                   %variance of the VAR constant (default=1e6)
hyperPriorsOptions.pos           =find(~hyperPars.isrw); %position of stationary variables
hyperPriorsOptions.MNalpha       =[];                    %lag decaying coeff of NIW prior (default=2)
hyperPriorsOptions.MNpsi         =false;                 %residual variance univariate AR(1) std (default=hyperprior)
hyperPriorsOptions.noc           =false;                 %sum of coefficients prior: NO default option
hyperPriorsOptions.sur           =false;                 %cointegration prior: NO default option
hyperPriorsOptions.Fcast         =false;                 %build forecasts: NO default option
hyperPriorsOptions.hz            =modelSpec.nHorizons;   %max forecast horizon: NO default option
hyperPriorsOptions.mcmc          =false;                 %run metropolis-hasting algorithm: NO default option
hyperPriorsOptions.Ndraws        =1200;                  %default=20k
hyperPriorsOptions.Ndrawsdiscard =200;                   %default=10k
hyperPriorsOptions.MCMCconst     =1;                     %default=1
hyperPriorsOptions.MCMCfcast     =false;                 %store forecast at each MCMC draw (default=true)
hyperPriorsOptions.MCMCstorecoeff=false;                 %store coefficients at each MCMC draw (default=true)
hyperPriorsOptions.initialValues =hyperPars;             %see above

hyperPriorsOptions.priorType = modelSpec.priorType;

%sampling from paramenters distribution
GibbsOptions.iterations          =1200;
GibbsOptions.burnin              =200;
GibbsOptions.jump                =1;

hyperPriorsOptions.GibbsOptions  =GibbsOptions;

%-------------------------------------------------------------------------%
%   ESTIMATE IMPULSE RESPONSE FUNCTIONS OVER DIFFERENT ROLLING SAMPLES    %
%-------------------------------------------------------------------------%


%variables in model
[~,selectN]                     =ismember(varList,VARdata.labels);

%data labels
macroVarNames                   =VARdata.description(selectN);
macroVarLabels                  =VARdata.labels(selectN);


%loop over samples
for t=1:modelSpec.nSubsamples
    

    %set estimation sample
    beginT                      =datenum(1965+t-1,1,1);          
    endT                        =datenum(1995+t-1,1,1);


    %set pre-sample for BLP
    beginPreT                   =datenum(1954+t-1,7,1); 
    endPreT                     =datenum(1964+t-1,10,1);
    
      
    %find relevant sample
    selectT                     =VARdata.dates>=beginT & VARdata.dates<=endT;
    selectPreT                  =VARdata.dates>=beginPreT & VARdata.dates<=endPreT;


    %use relevant observations

    macroVarData                =VARdata.data(selectT,selectN);
    macroVarDates               =VARdata.dates(selectT);


    macroVarPreData             =VARdata.data(selectPreT,selectN);
    macroVarPreDates            =VARdata.dates(selectPreT);


    %load data structure
    dataStructure.data          =macroVarData;
    dataStructure.dates         =macroVarDates;
    dataStructure.preSdata      =macroVarPreData;

    dataStructure.varname       =macroVarLabels;
    dataStructure.varLongName   =macroVarNames;

    
    %load shock variable and normalization in model structure 
    modelSpec.shockSize         =shockSize*double(ismember(dataStructure.varname,shockVar))';
    modelSpec.shockVar          =ismember(dataStructure.varname,shockVar);


    %load data in model specification
    modelSpec.dataStructure     =dataStructure;
    

    %  .     .     .     .     .     .     .     .     .     .     .     %
    %  .     .     .     .     .     .     .     .     .     .     .     %
    %  .     .     .     .     .     .     .     .     .     .     .     %

    %estimate & store IRFs over subsamples
    
    eval(['irfVAR', num2str(t),'=IRFbayesianNIW(modelSpec,hyperPriorsOptions);'])
    eval(['irfLP' , num2str(t),'=IRFlocalProj(modelSpec);'])
    eval(['irfBLP', num2str(t),'=IRFbayesianLocalProj(modelSpec,hyperPriorsOptions);'])

    %  .     .     .     .     .     .     .     .     .     .     .     %
    %  .     .     .     .     .     .     .     .     .     .     .     %
    %  .     .     .     .     .     .     .     .     .     .     .     %

    
end


%-SAVE OUTPUT-------------------------------------------------------------%

modelString='6519';

outFileName=['IRFs',modelString,'_',runIdentifier];
save([pwd '/MATfiles/' outFileName])



%%
%-PLOT IRFs---------------------------------------------------------------%

LineColors =[75 75 75;...
             0 50 100; %3
             0 50 100]./255;
          
BandColors =[.85 .85 .85;...
             .85 .85 .85;...
             [141 161 177]/255];

LineTypes  ={'-.',':'};

         
%plot labels
varname     =modelSpec.dataStructure.varLongName;
shockSize   =modelSpec.shockSize(modelSpec.shockVar)*100;
nHorizon    =modelSpec.nHorizons; 



%variables in plot
[~,selectN]  =ismember({'RGDP','RCON','HOUR','WAGE','FFR'},...
              modelSpec.dataStructure.varname);

%[~,selectN]=ismember(modelSpec.dataStructure.varname,modelSpec.dataStructure.varname);

plotColumns =numel(selectN);


figure; n=length(varname);

%-------BLP space----------------------------------------------------%
irfBLPareas=cell(n,1);
for j=selectN

    irfBLPall=nan(nHorizon+1,modelSpec.nSubsamples+1);
    
    %loop over subsamples
    for t=1:modelSpec.nSubsamples
        
        %load relevant IRF set
        eval(['irfBLPall(:,t)=irfBLP',num2str(t),'.irfs(:,j);'])
        
    end
    irfBLPareas{j}=[min(irfBLPall,[],2) max(irfBLPall,[],2)];
end        
          

pln=1;

%loop over variables
for j=selectN

    %  .     .     .     .     .     .     .     .     .     .     .     %
    %  .     .     .     .     .     .     .     .     .     .     .     %

    subplot(2,plotColumns,pln);

    hold on

    %BLP space
    p1=fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
        [irfBLPareas{j}(:,1)' fliplr(irfBLPareas{j}(:,2)')],...
        [141 161 177]/255,'EdgeColor','none');
    
    
    %loop over subsamples (VAR)
    for t=1:modelSpec.nSubsamples
        
        %load relevant IRF set
        eval(['irfVAR=irfVAR',num2str(t),';'])

        %VAR irfs
        p2=plot(0:nHorizon,irfVAR.irfs(:,j),LineTypes{1},'LineWidth',1,'color',LineColors(1,:));

        if t==modelSpec.nSubsamples

            %add zero line
            plot(0:nHorizon,zeros(size(1:nHorizon+1)),'r')

            %axis
            xlim([0 nHorizon]);     axis tight
            set(gca,'XTick',0:5:nHorizon,'XTickLabel',cellstr(num2str((0:5:nHorizon)')),'FontSize',9,'layer','top','Ygrid','on')
            title(varname{j},'FontSize',9,'FontWeight','normal')

        end        
    end
    
    
    if j==n    
        
        legend([p1 p2],{'BLP';'BVAR'},'FontSize',9,'Location','NorthEast','box','off')
    end

    %  .     .     .     .     .     .     .     .     .     .     .     %
    %  .     .     .     .     .     .     .     .     .     .     .     %
    
    subplot(2,plotColumns,pln+numel(selectN));

    hold on

    %BLP space
    p1=fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
        [irfBLPareas{j}(:,1)' fliplr(irfBLPareas{j}(:,2)')],...
        [141 161 177]/255,'EdgeColor','none');
    
    
    %loop over subsamples (LP)
    for t=1:modelSpec.nSubsamples
        
        %load relevant IRF set
        eval(['irfLP=irfLP',num2str(t),';'])

        %LP irfs
        p2=plot(0:nHorizon,irfLP.irfs(:,j),LineTypes{2},'LineWidth',1.1,'color',LineColors(2,:));

        if t==modelSpec.nSubsamples

            %add zero line
            plot(0:nHorizon,zeros(size(1:nHorizon+1)),'r')

            %axis
            xlim([0 nHorizon]);     axis tight
            set(gca,'XTick',0:5:nHorizon,'XTickLabel',cellstr(num2str((0:5:nHorizon)')),'FontSize',9,'layer','top','Ygrid','on')

        end        
    end
    
    
    if j==n    
        legend([p1 p2],{'BLP';'LP'},'FontSize',9,'Location','NorthEast','box','off')
    end
    
    
    pln=pln+1;
    
end

set(gcf,'PaperUnits','centimeters','PaperSize',[23 12]) %[x y]
set(gcf,'PaperPosition',[-2.5 0 27 12]) %[left bottom width height]
print(gcf,'-dpdf',[pwd '/CHARTS/' runIdentifier '_' modelString '.pdf']);            

