function res=DieboldMariano(fcstError1,fcstError2,horizon,squareLoss)

% Diebold&Mariano(1994) test for comparing predictive accuracy + DM test
% with modification for small sample as in Harvey,Leybourne&Newbold(1997)
%
% H0: equal predictive ability
% inputs: 2x [Tx1] forecast error series of competing models
%         forecast horizon (h); default truncation lag for sample
%         autocovariances set equal to (h-1)
%         1/0 for squareloss, if 0 loss=identity
%
% miranda (2013) silvia.mirandaagrippino@now-casting.com
% original code from R. Giacomini


T=length(fcstError1);

if T<10

    res.tstat   = nan;
    res.pval    = nan;
    res.tstat_ss= nan;
    res.pval_ss = nan;
    
    return
end



if squareLoss==1

    dt=fcstError1.^2-fcstError2.^2; %loss differential (squared loss fct)

else
    
    dt=fcstError1-fcstError2;
end

dbar=mean(dt); dt=dt-dbar;


%estimator for (+ve definite) asymptotic variance (Newey-West 1987)
Sigmat =var(dt,1);
lags   =horizon-1;

if lags>0
    for j=1:lags

        wj     =1-(j/(lags+1));
        Sigmat =Sigmat + wj*2*(dt(1+j:end)'*dt(1:end-j)/(T-j));
        
    end
end

tstat =dbar/sqrt(Sigmat/T); 
pval  =(1-normcdf(abs(tstat),0,1))*2;


%HLN small sample correction
tstat_ss =tstat/sqrt((T+1-2*horizon+horizon*(horizon-1)/T)/T); 
pval_ss  =(1-tcdf(abs(tstat_ss),T-1))*2;


%load results
res.tstat   = tstat;
res.pval    = pval;
res.tstat_ss= tstat_ss;
res.pval_ss = pval_ss;
