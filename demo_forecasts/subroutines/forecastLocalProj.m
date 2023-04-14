
function res=forecastLocalProj(modelSpec)

% produces LP-based univariate OOS forecasts for selected horizons
%
% inputs:
% modelSpec     =structure with model specification
%
% output: structure
% LPforecasts          =3D container of forecasted values
% LPforecastErrors     =3D contained of forecast errors
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
leveldata  =modelSpec.dataStructure.leveldata;
datafilter =modelSpec.dataStructure.tcodes;


%unpack forecast dates
startE     =find(dates==modelSpec.beginSample);
startF     =find(dates==modelSpec.beginForecast);
endF       =find(dates==modelSpec.endForecast);


%
n          =size(data,2);


%collectors
fcst_collect  =nan(endF-startF+max(Hstore),n,nH);
error_collect =nan(endF-startF+max(Hstore),n,nH);



%COMPUTE FORECASTS & FORECAST ERRORS

for j=1:n %for each variable
    
    for h=Hstore %for each horizon
        
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
                LPfcst =forecastLP(y,h,modelSpec);
            
                %revert to original data transform
                LPfcst =revert(LPfcst,h,leveldata(1:t,j),datafilter(j));                
                
                
                
                %forecasts
                fcst_collect(t+h-startF,j,Hstore==h) = LPfcst;
                
                %forecast errors
                error_collect(t+h-startF,j,Hstore==h)= leveldata(t+h,j) - fcst_collect(t+h-startF,j,Hstore==h);

            end
        end
    end    
end



%load results
res.forecasts   =fcst_collect;
res.errors      =error_collect;





%--------------------------------------------------------------------------
%-- child functions -------------------------------------------------------
%--------------------------------------------------------------------------


function res =forecastLP(y,h,modelSpec)

%univariate LP-based forecast: returns predicted values at selected
%horizon

%unpack model structure
nL    =modelSpec.nLPlags; %if lag selection this becomes the max lag allowed


%lag selection
if modelSpec.lagSelection
    
    nL =optimalLag(y,nL,h,'bic');
    
    %update model specification
    modelSpec.nLPlags =nL;
end


[T,n] =size(y); 


%build relevant projection set (shift obs backward to match horizon)
YhLag =nan(T-(nL+h),n*nL);

for j=h:nL+h-1

    YhLag(:,n*(j-h)+1:n*(j-h+1)) = y(nL+h-j+1:end-j,:);

end  

Yh        =y(nL+1+h:end,:); 
nT        =size(Yh,1);

YhprojSet =[ones(nT,1) YhLag];


%estimate local projection coefficients
projCoeffs =YhprojSet\Yh;


%forecast
res =[1 Yh(end:-1:end-nL+1,1)']*projCoeffs;





%--------------------------------------------------------------------------
function res =optimalLag(y,nL,h,criterion)

Ocrit =nan(nL,1);

[T,n] =size(y); 

l     =1;


while l <= nL 
        
    %build relevant projection set (shift obs backward to match horizon)
    YhLag =nan(T-(l+h),n*l);

    for j=h:l+h-1

        YhLag(:,n*(j-h)+1:n*(j-h+1)) = y(l+h-j+1:end-j,:);
    end  

    Yh        =y(l+1+h:end,:); 
    nT        =size(Yh,1);

    YhprojSet =[ones(nT,1) YhLag];


    switch criterion

        case 'aic'
            Ocrit(l) =aic(Yh,YhprojSet);

        case 'bic'
            Ocrit(l) =bic(Yh,YhprojSet);
    end
    
    if l>1 && find(Ocrit==min(Ocrit)) ~= sum(~isnan(Ocrit))
        
        break
    end
    l=l+1;
end

res =find(Ocrit==min(Ocrit));



%--------------------------------------------------------------------------
function res =revert(modelfcst,h,data,tcode)

%data at forecast origin: in levels
y     =data(end);
dy    =data(end)-data(end-1);

ly    =log(data(end));
dly   =log(data(end))-log(data(end-1));


%model forecast 
f_old =modelfcst;  %1 x 1


%revert transformation
switch tcode
   
    case 1 % level (i.e. no transformation): x(t)
        
        f_new =f_old;
        
        
    case 2 % first difference: x(t)-x(t-1)
        
        f_new =y + f_old; 
        
        
    case 3 % second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
        
        f_new =y + h*dy + f_old; 

        
    case 4 % natural log: ln(x)
        
        f_new =exp(f_old);
        
        
    case 5 % first difference of natural log: ln(x)-ln(x-1)
        
        f_new =exp(ly + f_old); 
        
        
    case 6 % second difference of natural log: (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
        
        f_new =exp(ly + h*dly + f_old); 
                
end


%load results
res = f_new;



