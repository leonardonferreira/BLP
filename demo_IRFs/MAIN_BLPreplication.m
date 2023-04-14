
% ----------------------------------------------------------------------- %
%                                                                         %
%                                                                         %
%          B A Y E S I A N   L O C A L   P R O J E C T I O N              %
%                                                                         %
%  Leonardo N. Ferreira, Silvia Miranda-Agrippino, Giovanni Ricco (2023)  %
%                                                                         %
%                                                                         %
%                  r e p l i c a t i o n   c o d e                        %
%                                                                         %
%                                                                         %
% Main script:                                                            %
% ------------                                                            %
%                                                                         %
% The following code estimates dynamic responses to a monetary policy     %
% shock on quarterly US data                                              % 
%                                                                         %
% TRANSMISSION:(options)                                                  %
%              - Local Projections                                        %
%              - Bayesian VAR                                             %
%              - Bayesian Local Projections with RW, VAR and DSGE priors  %
%                                                                         %
% IDENTIFICATION:(options)                                                %
%              - Cholesky                                                 %
%              - External Instruments                                     %
%                                                                         %
%                                                                         %
% ----------------------------------------------------------------------- %


clear
clc
load trueJPTirf;  
clearvars -except scaled_irf_y_c

start_time = clock;

addpath([pwd '/subroutines/']) %mac

%run identifier
%runIdentifier='JPTQ_BASELINE'; 
%runIdentifier='JPTQ_RWPRIOR';
runIdentifier='JPTQ_MODELBASED';



%-LOAD DATA---------------------------------------------------------------%

[tempData,tempText]      =xlsread('DataNew2019.xlsx','FRED_Data_Levels');
VARdata.data             =tempData(:,2:end);
VARdata.dates            =x2mdate(tempData(:,1));
VARdata.labels           =tempText(2,1:end);
VARdata.description      =tempText(1,1:end);



%-MODEL SPECIFICATION-----------------------------------------------------%


%choose VAR variables
varList                        =VARdata.labels;


%set estimation sample
beginT                         =datenum(1965,1,1);          
endT                           =datenum(2019,10,1);


%set pre-sample for BLP
beginPreT                      =datenum(1954,7,1); 
endPreT                        =datenum(1964,10,1);


%initialize BLP prior on presample
modelSpec.presample            =true;

%set BLP prior type
modelSpec.priorType            ='DSGE'; %BLP only, alternatives: 'VAR', 'RW', 'DSGE'; 
modelSpec.irfDSGE              =scaled_irf_y_c;


%set VAR specification, lags, IRFs horizon
modelSpec.nVARlags             =5; %number of lags in VAR (& prior)
modelSpec.nBLPlags             =5; %number of lags in BLP
modelSpec.nLPlags              =5; %number of lags in LP

modelSpec.nHorizons            =20;
modelSpec.bandsCoverage        =90;



%choose identification scheme
modelSpec.identification       ='CHOL'; %alternatives: 'PSVAR', 'CHOL';



%declare shock variable & shock size
shockVar                       ='FFR';       
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
%-------------------------------------------------------------------------%


%-LOAD DATA IN MODEL STRUCTURE--------------------------------------------%


%relevant sample
[~,selectN]                 =ismember(varList,VARdata.labels);
selectT                     =VARdata.dates>=beginT & VARdata.dates<=endT;
selectPreT                  =VARdata.dates>=beginPreT & VARdata.dates<=endPreT;


%use relevant observations
macroVarData                =VARdata.data(selectT,selectN);
macroVarDates               =VARdata.dates(selectT);


macroVarPreData             =VARdata.data(selectPreT,selectN);
macroVarPreDates            =VARdata.dates(selectPreT);


macroVarNames               =VARdata.description(selectN);
macroVarLabels              =VARdata.labels(selectN);


%load data structure
dataStructure.data          =macroVarData;
dataStructure.dates         =macroVarDates;
dataStructure.preSdata      =macroVarPreData;
dataStructure.preSdates     =macroVarPreDates;

dataStructure.varname       =macroVarLabels;
dataStructure.varLongName   =macroVarNames;


%load shock variable and normalization in model structure 
modelSpec.shockSize         =shockSize*double(ismember(dataStructure.varname,shockVar))';
modelSpec.shockVar          =ismember(dataStructure.varname,shockVar);


%load data in model specification
modelSpec.dataStructure     =dataStructure;

%-------------------------------------------------------------------------%
%                   ESTIMATE IMPULSE RESPONSE FUNCTIONS                   %
%-------------------------------------------------------------------------%


%estimate LP & compute IRFs
IRF_LP         =IRFlocalProj(modelSpec);

%estimate VAR & compute IRFs
IRF_BVAR       =IRFbayesianNIW(modelSpec,hyperPriorsOptions);

%estimate BLP & compute IRFs
IRF_BLP        =IRFbayesianLocalProj(modelSpec,hyperPriorsOptions);


%Including pre-sample in the VAR and in the classical LP
modelSpec.dataStructure.data  = [modelSpec.dataStructure.preSdata; modelSpec.dataStructure.data];
modelSpec.dataStructure.dates = [modelSpec.dataStructure.preSdates; modelSpec.dataStructure.dates];

IRF_LP2         =IRFlocalProj(modelSpec);
IRF_BVAR2       =IRFbayesianNIW(modelSpec,hyperPriorsOptions);


%-SAVE OUTPUT-------------------------------------------------------------%

modelString=[datestr(beginT,'YY'),datestr(endT,'YY')];

outFileName=['IRFs',modelString,'_',runIdentifier];
save([pwd '/MATfiles/' outFileName])

%-------------------------------------------------------------------------&



disp( ['Estimation takes '  num2str( etime( clock, start_time) ) ' seconds' ] );



%%
%-PLOT IRFs---------------------------------------------------------------%

LineColors =[75 75 75;...
             75 75 75; %3
             0 50 100]./255;
          
BandColors =[.85 .85 .85;...
             .85 .85 .85;...
             [141 161 177]/255];

LineTypes  ={'-.','--','x'};
         
%plot labels
varname     =modelSpec.dataStructure.varLongName;
shockSize   =modelSpec.shockSize(modelSpec.shockVar)*100;
nHorizon    =modelSpec.nHorizons; 



%variables in plot
% [~,selectN]  =ismember({'RGDP','RCON','HOUR','WAGE','FFR'},...
%               modelSpec.dataStructure.varname);

[~,selectN]=ismember(modelSpec.dataStructure.varname,modelSpec.dataStructure.varname);


plotColumns =numel(selectN);


%-------BLP VS VAR AND LP-------------------------------------------------%
figure; pln=1; n=length(varname);
    
for j=selectN
    
    %  .     .     .     .     .     .     .     .     .     .     .     %
    %  .     .     .     .     .     .     .     .     .     .     .     %

    %BLP vs VAR
    subplot(2,plotColumns,pln);
                    
    hold on
        
    %bands
%     fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
%         [IRF_BVAR.irfs_l(1:nHorizon+1,j)' fliplr(IRF_BVAR.irfs_u(1:nHorizon+1,j)')],...
%         BandColors(2,:),'EdgeColor','none');

    fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
        [IRF_BLP.irfs_l(1:nHorizon+1,j)' fliplr(IRF_BLP.irfs_u(1:nHorizon+1,j)')],...
        BandColors(3,:),'EdgeColor','none');
    
    fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
        [IRF_BVAR.irfs_l(1:nHorizon+1,j)' fliplr(IRF_BVAR.irfs_u(1:nHorizon+1,j)')],...
        BandColors(2,:),'EdgeColor','none');

    %irfs
    p1=plot(0:nHorizon,IRF_BVAR.irfs(1:nHorizon+1,j),'--','LineWidth',1.2,'color',LineColors(1,:));
        
    p2=plot(0:nHorizon,IRF_BLP.irfs(1:nHorizon+1,j),'x','LineWidth',1.2,'color',LineColors(3,:));

    
    %zero line
    plot(0:nHorizon,zeros(size(1:nHorizon+1)),'r','LineWidth',.7)
    
    
    %axis
    xlim([0 nHorizon]); axis tight
    yticks=get(gca,'Ylim'); yticks=round(linspace(yticks(1), yticks(end), 5),1);
    set(gca,'YTick',yticks,'XTick',0:5:nHorizon,'XTickLabel',cellstr(num2str((0:5:nHorizon)')),'FontSize',11,'Layer','top','Ygrid','on')
    title(varname{j},'FontSize',13,'FontWeight','Normal')

    
    
    if j==max(selectN)
            
        lh=legend([p1 p2],{['BVAR(' num2str(modelSpec.nVARlags) ')'],...
                             ['BLP(' num2str(modelSpec.nBLPlags) ')']},'FontSize',11,'Location','NorthEast','box','off');

    
    end         
           
    %BLP vs LP
    subplot(2,plotColumns,pln+numel(selectN));
                    
    hold on
        
    %bands
    fill([1:1:nHorizon, fliplr(1:1:nHorizon)],...
        [IRF_LP.irfs_l(2:nHorizon+1,j)' fliplr(IRF_LP.irfs_u(2:nHorizon+1,j)')],...
        BandColors(1,:),'EdgeColor','none');
    
    fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
        [IRF_BLP.irfs_l(1:nHorizon+1,j)' fliplr(IRF_BLP.irfs_u(1:nHorizon+1,j)')],...
        BandColors(3,:),'EdgeColor','none');
    
    %irfs
    p3=plot(0:nHorizon,IRF_LP.irfs(1:nHorizon+1,j),'-.','LineWidth',1.2,'color',LineColors(2,:));
    
    plot(0:nHorizon,IRF_BLP.irfs(1:nHorizon+1,j),'x','LineWidth',1.2,'color',LineColors(3,:));
    
    
    %zero line
    plot(0:nHorizon,zeros(size(1:nHorizon+1)),'r','LineWidth',.7)
    
    
    %axis
    xlim([0 nHorizon]); axis tight
    yticks=get(gca,'Ylim'); yticks=round(linspace(yticks(1), yticks(end), 5),1);
    set(gca,'YTick',yticks,'XTick',0:5:nHorizon,'XTickLabel',cellstr(num2str((0:5:nHorizon)')),'FontSize',11,'Layer','top','Ygrid','on')
    title(varname{j},'FontSize',13,'FontWeight','Normal')
        
    if j==max(selectN)
        
        lh=legend([p3 p2],{['LP(' num2str(modelSpec.nLPlags) ')'],...
            ['BLP(' num2str(modelSpec.nBLPlags) ')']},'FontSize',11,'Location','NorthEast','box','off');
        
    end
    
    
    pln=pln+1;
end


set(gcf,'PaperUnits','centimeters','PaperSize',[29 18]) %[x y]
set(gcf,'PaperPosition',[-3.5 0 35 18]) %[left bottom width height]


print(gcf,'-dpdf',[pwd '/CHARTS/' runIdentifier '_' modelString '_Bands.pdf']);
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------BLP VS VAR AND LP 2 ----------------------------------------------%
%VAR and LP use training sample
figure; pln=1; n=length(varname);
    
for j=selectN
    
    %  .     .     .     .     .     .     .     .     .     .     .     %
    %  .     .     .     .     .     .     .     .     .     .     .     %

    %BLP vs VAR
    subplot(2,plotColumns,pln);
                    
    hold on
        
    %bands
%     fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
%         [IRF_BVAR.irfs_l(1:nHorizon+1,j)' fliplr(IRF_BVAR.irfs_u(1:nHorizon+1,j)')],...
%         BandColors(2,:),'EdgeColor','none');

    fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
        [IRF_BLP.irfs_l(1:nHorizon+1,j)' fliplr(IRF_BLP.irfs_u(1:nHorizon+1,j)')],...
        BandColors(3,:),'EdgeColor','none');
    
    fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
        [IRF_BVAR2.irfs_l(1:nHorizon+1,j)' fliplr(IRF_BVAR2.irfs_u(1:nHorizon+1,j)')],...
        BandColors(2,:),'EdgeColor','none');

    %irfs
    p1=plot(0:nHorizon,IRF_BVAR2.irfs(1:nHorizon+1,j),'--','LineWidth',1.2,'color',LineColors(1,:));
        
    p2=plot(0:nHorizon,IRF_BLP.irfs(1:nHorizon+1,j),'x','LineWidth',1.2,'color',LineColors(3,:));

    
    %zero line
    plot(0:nHorizon,zeros(size(1:nHorizon+1)),'r','LineWidth',.7)
    
    
    %axis
    xlim([0 nHorizon]); axis tight
    yticks=get(gca,'Ylim'); yticks=round(linspace(yticks(1), yticks(end), 5),1);
    set(gca,'YTick',yticks,'XTick',0:5:nHorizon,'XTickLabel',cellstr(num2str((0:5:nHorizon)')),'FontSize',10,'Layer','top','Ygrid','on')
    title(varname{j},'FontSize',11,'FontWeight','normal')

    
    
    if j==max(selectN)
            
        lh=legend([p1 p2],{['BVAR(' num2str(modelSpec.nVARlags) ')'],...
                             ['BLP(' num2str(modelSpec.nBLPlags) ')']},'FontSize',11,'Location','Northeast','box','off');
    
    end         
    
        
    %BLP vs LP
    subplot(2,plotColumns,pln+numel(selectN));
                    
    hold on
        
    %bands
    fill([1:1:nHorizon, fliplr(1:1:nHorizon)],...
        [IRF_LP2.irfs_l(2:nHorizon+1,j)' fliplr(IRF_LP2.irfs_u(2:nHorizon+1,j)')],...
        BandColors(1,:),'EdgeColor','none');
    
    fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
        [IRF_BLP.irfs_l(1:nHorizon+1,j)' fliplr(IRF_BLP.irfs_u(1:nHorizon+1,j)')],...
        BandColors(3,:),'EdgeColor','none');
    
    %irfs
    p3=plot(0:nHorizon,IRF_LP2.irfs(1:nHorizon+1,j),'-.','LineWidth',1.2,'color',LineColors(2,:));
    
    plot(0:nHorizon,IRF_BLP.irfs(1:nHorizon+1,j),'x','LineWidth',1.2,'color',LineColors(3,:));
    
    
    %zero line
    plot(0:nHorizon,zeros(size(1:nHorizon+1)),'r','LineWidth',.7)
    
    
    %axis
    xlim([0 nHorizon]); axis tight
    yticks=get(gca,'Ylim'); yticks=round(linspace(yticks(1), yticks(end), 5),1);
    set(gca,'YTick',yticks,'XTick',0:5:nHorizon,'XTickLabel',cellstr(num2str((0:5:nHorizon)')),'FontSize',10,'Layer','top','Ygrid','on')
    title(varname{j},'FontSize',11,'FontWeight','normal')
    
    
    if j==max(selectN);
        
        lh=legend([p3 p2],{['LP(' num2str(modelSpec.nLPlags) ')'],...
            ['BLP(' num2str(modelSpec.nBLPlags) ')']},'FontSize',11,'Location','Northeast','box','off');
        
    end
    
    
    pln=pln+1;
end

set(gcf,'PaperUnits','centimeters','PaperSize',[29 18]) %[x y]
set(gcf,'PaperPosition',[-3.5 0 35 18]) %[left bottom width height]

print(gcf,'-dpdf',[pwd '/CHARTS/' runIdentifier '_' modelString '_Bands_FS.pdf']);
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

%%

%-------BLP OPTIMAL TIGHTNESS-----------------------------------------%
figure;

hold on

%h=1 (VAR)
plot([IRF_BLP.optimalLambda(1); nan(nHorizon-1,1)],'s','LineWidth',1.2,'color',[141 161 177]/255,'MarkerSize',8,'MarkerFaceColor',[141 161 177]/255);
%h>1 (BLP)
plot([nan; IRF_BLP.optimalLambda(2:nHorizon)]   ,'s','LineWidth',1.2,'color',LineColors(3,:),'MarkerSize',8,'MarkerFaceColor',LineColors(3,:));

hold off

%axis
set(gca,'FontSize',10,'Ygrid','on'); xlim([0 nHorizon+1]); %grid on
xlabel('horizon','FontSize',11)


set(gcf,'PaperUnits','centimeters','PaperSize',[20 9]) %[x y]
set(gcf,'PaperPosition',[-1 0 22 9]) %[left bottom width height]
print(gcf,'-dpdf',[pwd '/CHARTS/' runIdentifier '_' modelString '_Shrinkage.pdf']);
