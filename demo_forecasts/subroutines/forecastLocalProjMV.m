
function res=forecastLocalProjMV(modelSpec)

% produces LP-based multivariate OOS forecasts for selected horizons
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

for h=Hstore %for each horizon

    for t=startF:endF %for each forecast origin

        %find start date for estimation
        switch fcstType

            case 'recursive'
                beginS =startE;      %sample alway starts from same obs

            case 'rolling'
                beginS =t-startF+1;  %length fixed to first available sample

        end


        y=data(beginS:t,:);


        %produce forecasts
        LPfcst =forecastLP(y,h,modelSpec);


        %forecasts
        fcst_collect(t+h-startF,:,Hstore==h) = LPfcst;
        
        %forecast errors
        error_collect(t+h-startF,:,Hstore==h)= 1/h*(data(t+h,:) - LPfcst);

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



%last set of obs
YTlag=y(end:-1:end-nL+1,:)';
YTlag=YTlag(:);


%forecast
res =[1 YTlag']*projCoeffs;



