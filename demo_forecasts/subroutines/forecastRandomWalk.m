
function res=forecastRandomWalk(modelSpec)

% produces RW-based univariate OOS forecasts for selected horizons
%
% inputs:
% modelSpec     =structure with model specification
%
% output: structure
% RWforecasts          =3D container of forecasted values
% RWforecastErrors     =3D contained of forecast errors
%
%
% miranda 2019 silvia.miranda-agrippino@northwestern.edu

%--------------------------------------------------------------------------


%unpack model structure
Hstore     =modelSpec.nHorizons;
nH         =numel(Hstore);            %forecast horizons

fcstType   =modelSpec.forecastType;   %recursive/rolling


%unpack data
data       =modelSpec.dataStructure.data;
dates      =modelSpec.dataStructure.dates;
% datafilter =modelSpec.dataStructure.tcodes;


%unpack forecast dates
startE     =find(dates>=modelSpec.beginSample,1,'first');
startF     =find(dates>=modelSpec.beginForecast,1,'first');
endF       =find(dates<=modelSpec.endForecast,1,'last');


%
n          =size(data,2);


%collectors
fcst_collect  =nan(endF-startF+max(Hstore),n,nH);
error_collect =nan(endF-startF+max(Hstore),n,nH);



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


        if modelSpec.dataStructure.ragged     %data series have different start date            
            if numel(y)<=120 || any(isnan(y)) %at least 10 years of data & no nans
                
                continue;
            end
        end
            
        %produce forecasts
        RWfcst =forecastRW(y,modelSpec);


%         %revert to original data transform
%         RWfcst =revert(RWfcst,leveldata(1:t,j),datafilter(j));



        for h=Hstore %for each horizon
           

            %forecasts
            fcst_collect(t+h-startF,j,Hstore==h) = RWfcst(:,h);

            %forecast errors
            error_collect(t+h-startF,j,Hstore==h)= 1/h*(data(t+h,j) - RWfcst(:,h));

        end    
    end
end



%load results
res.forecasts   =fcst_collect;
res.errors      =error_collect;





%--------------------------------------------------------------------------
%-- child functions -------------------------------------------------------
%--------------------------------------------------------------------------


function res =forecastRW(y,modelSpec)

%univariate RW-based forecast: returns predicted values at selected
%horizon


%unpack model structure
nH   =max(modelSpec.nHorizons);
n    =size(y,2); 


%store forecasts
fcst =nan(n,nH);

for h=1:nH

    fcst(:,h) =y(end,:) + (y(end,:)-y(end-h,:)); %RW=constant growth
end

%load results
res=fcst;


%--------------------------------------------------------------------------
function res =revert(modelfcst,data,tcode)

%data at forecast origin: in levels
y     =data(end);
dy    =data(end)-data(end-1);

ly    =log(data(end));
dly   =log(data(end))-log(data(end-1));


%model forecast 
f_old =modelfcst;  %1 x nH

%max forecast horizon
nH    =numel(f_old);



%revert transformation
%random walk is treated as a direct forecast
switch tcode
   
    case 1 % level (i.e. no transformation): x(t)
        
        f_new =f_old;
        
        
    case 2 % first difference: x(t)-x(t-1)
        
        f_new =nan(1,nH);        
        for h=1:nH
           
            f_new(h) =y + f_old(h); 
        end        
        
        
    case 3 % second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
        
        f_new =nan(1,nH);        
        for h=1:nH
           
            f_new(h) =y + h*dy + f_old(h); 
        end        
        
        
    case 4 % natural log: ln(x)
        
        f_new =exp(f_old);
        
        
    case 5 % first difference of natural log: ln(x)-ln(x-1)
        
        f_new =nan(1,nH);        
        for h=1:nH
           
            f_new(h) =ly + f_old(h); 
        end     
        f_new =exp(f_new);
        
        
    case 6 % second difference of natural log: (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
        
        f_new =nan(1,nH);        
        for h=1:nH
           
            f_new(h) =ly + h*dly + f_old(h); 
        end
        f_new =exp(f_new);
                
end


%load results
res = f_new;


