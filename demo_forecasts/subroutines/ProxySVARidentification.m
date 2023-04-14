
% PROXY SVAR with instrumental variable identification of IRFs
% built on Mertens and Ravn (2012)
%
% miranda 2014 smirandaagrippino@london.edu
%
%
%
% Modified to allow identification of multiple shocks (m>1) with multiple
% (m) instruments which may be cross-correlated. 
% Assumes the variables enter the VAR in the desired Cholesky ordering
%
% miranda 2016 silvia.miranda-agrippino@bankofengland.co.uk


function res = ProxySVARidentification(VARinnovations,policyVarPos,proxyVar)

% inputs:
% VARinnovations = [nxT] matrix of reduced form innovations
% policyVarPos   = position of the policy indicators within system
% iVar           = [Txm] vector of instruments for identification of policy shocks
%
% output (structure)
% B              = [Nxm] vector of impact responses of VAR innovations to structural policy shocks
% psi            = estimated correlation between instument and structural policy shock
% L              = reliability of instrument
% e              = [1xT] vector of realized policy shock
% Fstat          = for the regression of the policy innovation on the instrument.



[n,T] =size(VARinnovations); %system size
m     =size(proxyVar,2);     %number of instruments


if numel(policyVarPos) ~= m
   
    error('the number of instruments must equal the number of shocks to be identified \n')
end


%place policy indicators first
iP  =policyVarPos;      niP =1:n; niP(iP)=[];

VARu    =[VARinnovations(iP,:);VARinnovations(niP,:)]; 
VARu    =VARu'; %Txn


% -coefficients of regression on instrument - - - - - - - - - - - - - - - %

betaIV =kron( eye(n), [ones(T,1) proxyVar] )\VARu(:);
betaIV =reshape(betaIV,length(betaIV)/n,n)'; %nx(m+1)

beta_11  =betaIV(1:m,2:m+1);
beta_21  =betaIV(m+1:end,2:m+1);

%ratio of regression coeffs 
B21B11   =beta_21/beta_11; %beta_{21}*beta_{11}^{-1}


% -identification (notation ~ as in Mertens and Ravn) - - - - - - - - - - %

SigmaU =cov(VARu);

Zeta   =B21B11*SigmaU(1:m,1:m)*B21B11' -...
        ( SigmaU(m+1:end,1:m)*B21B11' + B21B11*SigmaU(m+1:end,1:m)' ) +...
        SigmaU(m+1:end,m+1:end);

B12B12 =( SigmaU(m+1:end,1:m) - B21B11*SigmaU(1:m,1:m) )'/Zeta*...
        ( SigmaU(m+1:end,1:m) - B21B11*SigmaU(1:m,1:m) ); %B12B12'


B11B11 =SigmaU(1:m,1:m) - B12B12; %B11B11'


if m ==1 %one instrument for one structural shock
    
    B11 =sqrt( B11B11 );   %\beta_{11}
    
    B   =B11.*[1; B21B11]; %first column of B (u_t = B*e_t)
       
    
else %m instruments for m shocks
    
    B22B22 =SigmaU(m+1:end,m+1:end) +...
            B21B11*( B12B12 - SigmaU(1:m,1:m) )*B21B11'; %B22B22'
    
    B12B22 =( B12B12*B21B11' + (SigmaU(m+1:end,1:m) - B21B11*SigmaU(1:m,1:m))' )/B22B22; %B12(B22^-1)
    
    B11S1  =eye(m) - B12B22*B21B11; %B11(S1^-1)
    
    B21S1  =B21B11/B11S1; %B21(S1^-1)
    
    S1S1   =B11S1*B11B11*B11S1';
    
    %variables enter the VAR in relevant Cholesky ordering
    S1     =chol(S1S1)';
    
    B      =[ inv(B11S1); B21S1 ]*S1; %first m columns of B (u_t = B*e_t)

end


%F stat (regression on instruments of relevant innovations)
tempX  =[ones(T,1) proxyVar]; 
tempU  =VARu(:,1:m) - tempX*betaIV(1:m,:)'; 

tempY  =tempX*betaIV(1:m,:)' - repmat(mean(VARu(:,1:m)),T,1); k = length(betaIV(1:m,:)')-1;
F_Stat =( (tempY'*tempY)/k )/( (tempU'*tempU)/(T-k-1) );



%realized shock sequences (Montiel-Olea, Stock and Watson)
tempX =[ones(T,1) VARu]; 
e     =tempX*(tempX\proxyVar);
e     =bsxfun(@rdivide,e,std(e));  %unit variance shock series


%reliability
if m==1
    
    D       =proxyVar~=0; %proportion of uncensored data

    SigmaMU =cov([proxyVar(D),VARu(D,:)],1); 

    %relevance
    Phi   =SigmaMU(1,2:end)/B';
    Gamma =(sum(D)/T)\Phi;

    
    eSquare =e.^2; zSquare =(proxyVar-Gamma*e).^2;

    Lambda  =( Gamma^2*sum(eSquare(D)) + sum(zSquare(D)) )\ Gamma^2*sum(eSquare(D));
    %v=proxyVar-Gamma*e; SigmaV=cov(v(D),1); Lambda=(Gamma^2+SigmaV)\Gamma^2; %alternative
    
else
    
    SigmaMU =cov([proxyVar,VARu(:,1:m)],1); 

    SigmaMM =cov(proxyVar,1);
    
    D =sum(proxyVar,2)~=0; %proportion of uncensored data    
    Lambda  =kron(T/sum(D),eye(m))/SigmaMM*SigmaMU(1:m,m+1:end)/B11B11*SigmaMU(m+1:end,1:m); 
    
    %relevance
    Gamma   =sqrt(sort(eig(Lambda)));

end


%load output structure
res.B        =NaN(n,m);      %contemporaneous transmission coefficients: Bzero

res.B(iP,:)  =B(1:m,1:m); 
res.B(niP,:) =B(m+1:end,:);
    
    
res.Gamma    =Gamma;         %estimated correlation between shock and instrument
res.L        =Lambda;        %reliability of instrument
res.e        =e;             %realized shocks series
res.fstat    =diag(F_Stat);  %F statistic of regression on instrument



