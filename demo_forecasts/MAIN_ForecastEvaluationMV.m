
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
% Forecast Evaluation                                                     %
% -------------------                                                     %
%                                                                         %                                                                        %
% MODELS:                                                                 %
%              - Random Walk                                              %
%              - Local Projections                                        %
%              - Bayesian VAR                                             %
%              - Bayesian Local Projections with RW, VAR and DSGE priors  %
%                                                                         %
% ----------------------------------------------------------------------- %

clear
clc


addpath([pwd '/subroutines/']) %mac
addpath([pwd '/MATfiles/']) %mac


%load forecasts
load FCST_1990-2017_MULTIVARIATE_RW_recursive_paper



% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %
%basics
n       =size(modelSpec.dataStructure.data,2);    %number of variables
nH      =numel(modelSpec.nHorizons);              %forecast horizons
nL      =modelSpec.nLPlags;                       %lags in LP

nM      =4;                                       %number of models

labels  =modelSpec.dataStructure.varname;         %variable names


startF  =find(modelSpec.dataStructure.dates>=modelSpec.beginForecast,1,'first');
endF    =find(modelSpec.dataStructure.dates<=modelSpec.endForecast,1,'last');



% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %
%1. POINT FORECASTS (ALL MODELS)

%RMSFE


%LP
RMSFE_RW   =squeeze(sqrt(nanmean(FCST_RW.errors.^2,1)));  

%LP
RMSFE_LP   =squeeze(sqrt(nanmean(FCST_LP.errors.^2,1)));  

%BVAR
RMSFE_BVAR =squeeze(sqrt(nanmean(FCST_VAR.errors.^2,1))); 

%BLP
RMSFE_BLP  =squeeze(sqrt(nanmean(FCST_BLP.errors.^2,1))); 


%accuracy relative to RW

relativeRMSFE_LP =RMSFE_LP./RMSFE_RW;
relativeRMSFE_BVAR =RMSFE_BVAR./RMSFE_RW;
relativeRMSFE_BLP =RMSFE_BLP./RMSFE_RW;


% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %
%                        TABLE 1: Average RMSFE
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %
TABLE1 =nan(n,nH*nM);

for h=1:nH
   
    TABLE1(:,nM*(h-1)+1:nM*h) =[RMSFE_RW(:,h) RMSFE_LP(:,h) RMSFE_BVAR(:,h) RMSFE_BLP(:,h)];
end
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %


%test of equal predictive accuracy (DM, for each variable, against RW)
%H0: equal predictive ability
DMtest =nan(n,nH,3);                %on third dimension: LP, VAR, BLP in this order

for j=1:n
   for h=1:nH
       
       %LP
       DM_LP         =DieboldMariano( FCST_LP.errors(~isnan(FCST_LP.errors(:,j,h)),j,h),...
                      FCST_RW.errors(~isnan(FCST_RW.errors(:,j,h)),j,h),...
                      modelSpec.nHorizons(h),1 );
           
       DMtest(j,h,1) =DM_LP.pval_ss;
       
       
              
       %VAR
       DM_VAR         =DieboldMariano( FCST_VAR.errors(~isnan(FCST_VAR.errors(:,j,h)),j,h),...
                      FCST_RW.errors(~isnan(FCST_RW.errors(:,j,h)),j,h),...
                      modelSpec.nHorizons(h),1 );
           
       DMtest(j,h,2) =DM_VAR.pval_ss;
       
       
       
       %BLP
       DM_BLP        =DieboldMariano( FCST_BLP.errors(~isnan(FCST_BLP.errors(:,j,h)),j,h),...
                      FCST_RW.errors(~isnan(FCST_RW.errors(:,j,h)),j,h),...
                      modelSpec.nHorizons(h),1 );
                  
       DMtest(j,h,3) =DM_BLP.pval_ss;
   end    
end



% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %
%                        TABLE 1b: Average Relative RMSFE
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %

TABLE1b =nan(2*n,nH*(nM-1));

for j=1:n
    
    TABLE1b(2*(j-1)+1:2*(j),:) = [relativeRMSFE_LP(j,:) relativeRMSFE_BVAR(j,:) relativeRMSFE_BLP(j,:);
                                     DMtest(j,:,1)      DMtest(j,:,2)  DMtest(j,:,3)    ];
end


% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %


%%

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %
%2. DENSITY FORECASTS (BVAR) 

nT       =size(FCST_VAR.errors,1);
data     =modelSpec.dataStructure.data;         %already transformed
horizon  =modelSpec.nHorizons;


%collectors
PITs_VAR =nan(nT,n,nH);
LogS_VAR =nan(nT,n,nH);


for j=1:n
    for h=1:nH
        for t=startF:endF %for each forecast origin
            
            
           
            
            VAR_density         =squeeze(FCST_VAR.densityf(t+horizon(h)-startF,j,h,:)); % T x n x H x S
            
            if ~all(isnan(VAR_density))
                              
                PITs_VAR(t+horizon(h)-startF,j,h) =ksdensity(VAR_density,data(t+horizon(h),j),...
                    'function','cdf','kernel','epanechnikov');
                
                LogS_VAR(t+horizon(h)-startF,j,h) =log(ksdensity(VAR_density,data(t+horizon(h),j),...
                    'function','pdf','kernel','epanechnikov'));
            end
                       
             
        end
    end
end

LogS_VAR(isinf(LogS_VAR)) =nan;


%% .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
%3. DENSITY FORECASTS with LP & BLP

LP_FEVariance =nan(n,nH); % variance of LP forecast errors (corrected for autocorrelation)
BLP_FEVariance =nan(n,nH); % variance of BLP forecast errors (corrected for autocorrelation)

for h=1:nH   
        
    %LP
    %forecast errors
    u           =FCST_LP.errors(:,:,h);   
    
    %remove nans (container is larger than number of forecasts)
    keeprows    =~any(isnan(u'));
    u           =u(keeprows,:);
    u           =bsxfun(@minus,u,mean(u)); 
    
    tau         =size(u,1); %length of forecasts sequence
       
    SigmaU      =u'*u; %HAC error (T*)covariance estimator

    nwLags      =h+1;
    
    nwWeights   =(nwLags+1-(1:nwLags))./(nwLags+1);
    for j=1:nwLags

        Gammaj  =(u(j+1:tau,:)'*u(1:tau-j,:));
        SigmaU  =SigmaU+nwWeights(j)*(Gammaj+Gammaj');

    end   
   
    %store
    LP_FEVariance(:,h)=diag(SigmaU)/(tau-nL);
 
    
    %BLP
    %forecast errors
    u           =FCST_BLP.errors(:,:,h);   
    
    %remove nans (container is larger than number of forecasts)
    keeprows    =~any(isnan(u'));
    u           =u(keeprows,:);
    u           =bsxfun(@minus,u,mean(u)); 
    
    tau         =size(u,1); %length of forecasts sequence
       
    SigmaU2      =u'*u; %HAC error (T*)covariance estimator

    nwLags      =h+1;
    
    nwWeights   =(nwLags+1-(1:nwLags))./(nwLags+1);
    for j=1:nwLags

        Gammaj2  =(u(j+1:tau,:)'*u(1:tau-j,:));
        SigmaU2 =SigmaU2+nwWeights(j)*(Gammaj2+Gammaj2');

    end   
   
    %store
    BLP_FEVariance(:,h)=diag(SigmaU2)/(tau-nL);
    
    
end



%Log predictive scores for LP & BLP (assumes normality)
LogS_LP =nan(nT,n,nH);
LogS_BLP =nan(nT,n,nH);

%calculate log scores
for j=1:n
    for h=1:nH
        for t=1:nT
        
            LogS_LP(t,j,h) =log(normpdf(h*FCST_LP.errors(t,j,h)+FCST_LP.forecasts(t,j,h),...
                                        FCST_LP.forecasts(t,j,h),sqrt(LP_FEVariance(j,h))));  
                                    
            LogS_BLP(t,j,h)=log(normpdf(h*FCST_BLP.errors(t,j,h)+FCST_BLP.pointf(t,j,h),...
                                        FCST_BLP.pointf(t,j,h),sqrt(BLP_FEVariance(j,h))));  
        end
    end
end
          

%%

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %
%                  TABLE 2: Average Difference in LS
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %

TABLE2 =nan(2*n,2*nH);

for j=1:n
        

    % difference in LS (per variable across horizons)
    diff_1 =squeeze(LogS_BLP(:,j,:) - LogS_LP(:,j,:));
    diff_2 =squeeze(LogS_BLP(:,j,:) - LogS_VAR(:,j,:));


    %load in table
    TABLE2(2*(j-1)+1:2*(j),:) = [nanmean(diff_1) nanmean(diff_2); nanstd(diff_1) nanstd(diff_2)];

end


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

%bands coverage
bandSize=68;

uBound  =bandSize+(100-bandSize)/2;   uBound=uBound/100;
lBound  =(100-bandSize+1)/2;          lBound=lBound/100;
sLevel=abs(norminv((1-bandSize/100)/2,0,1));

nD      =hyperPriorsOptions.GibbsOptions.iterations-hyperPriorsOptions.GibbsOptions.burnin;

horizon =modelSpec.nHorizons;

%plots (density)
figure;         pln=1; 

for h=1:H

    for j=1:n
        
    
        VAR_density         =squeeze(FCST_VAR.densityf(:,j,h,:)./100);
        BLP_upper           = (squeeze(FCST_BLP.pointf(:,j,h)) + BLP_FEVariance(j,h).^0.5*sLevel)/100;
        BLP_lower           = (squeeze(FCST_BLP.pointf(:,j,h)) - BLP_FEVariance(j,h).^0.5*sLevel)/100;
        
    
        subplot(H,n,pln)
        
        hold on

        %bands 
        p1=fill([dates(horizon(h):end-horizon(H)+horizon(h)-1)', fliplr(dates(horizon(h):end-horizon(H)+horizon(h)-1)')],...
            [VAR_density(horizon(h):end-horizon(H)+horizon(h)-1,round(uBound*nD))' fliplr(VAR_density(horizon(h):end-horizon(H)+horizon(h)-1,round(lBound*nD))')],...
            [186 196 196]/255,'EdgeColor','none');
        
        %bands 
        p2=fill([dates(horizon(h):end-horizon(H)+horizon(h)-1)', fliplr(dates(horizon(h):end-horizon(H)+horizon(h)-1)')],...
            [BLP_upper(horizon(h):end-horizon(H)+horizon(h)-1)' fliplr(BLP_lower(horizon(h):end-horizon(H)+horizon(h)-1)')],...
            [141 161 177]/255,'EdgeColor','none');

        p3=plot(dates,squeeze(FCST_VAR.pointf(:,j,h))./100,'-.','color',[141 161 177]/255);
        p4=plot(dates,squeeze(FCST_BLP.pointf(:,j,h))./100,'-','color',[0 50 100]/255);
        
        p5=plot(dates,modelSpec.dataStructure.data(startF:endF+horizon(H)-1,j)./100,'-','color',Ocolor);
                
        
        set(gca,'FontSize',7,'Xtick',dates(1:36:end),'box','on','layer','top','Ygrid','on'); 
        dateaxis('x',10); datetick('keepticks'); xlim([dates(1) dates(end)]);
        title(labels{j},'Fontweight','Normal','Fontsize',8)
        if j==1 && h==1
            legend([p5 p1 p2],{'data','VAR','BLP'},'FontSize',8,'Location','SouthEast','box','off')
        end
        
        pln=pln+1;
    end
end

set(gcf,'PaperUnits','centimeters','PaperSize',[29 18]) %[x y]
set(gcf,'PaperPosition',[-3.5 0 35 18]) %[left bottom width height]
print(gcf,'-dpdf',[pwd '/CHARTS/DensityForecasts_' modelString '.pdf']);            



