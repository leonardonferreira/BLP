function omegahat = NeweyWest(Z,nlags)

% Returns the Newey-West estimator of the asymptotic variance matrix
% INPUTS: Z, a nxk matrix with rows the vector zt'
%         nlags, the number of lags
%
% OUTPUTS: omegahat, the Newey-West estimator of the covariance matrix


[n,k] = size(Z);

% de-mean the variables
Z = Z - ones(size(Z,1),1)*mean(Z);

gamma = -999*ones(nlags,k);
samplevar = Z'*Z/n; % sample variance
omegahat = samplevar;
if nlags > 0
   % sample autocovariances
   for ii = 1:nlags
      Zlag = [zeros(ii,k);Z(1:n-ii,:)];
      gamma = (Z'*Zlag +Zlag'*Z)/n;
      weights = 1 - (ii/(nlags+1));
      omegahat = omegahat + weights*gamma;
   end
end



   
   










