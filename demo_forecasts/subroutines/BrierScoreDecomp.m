function [BrierScore] = BrierScoreDecomp(probForecast, eventOccured,nBins)

nRows               = size(probForecast, 1);
% Initialise bins and occ
binUnitLocation = 0: 1 / nBins : 1;
Bin(1:nBins)    = struct('bin', 0);
Occ(1:nBins)    = struct('occ', 0);

% Need adjusting for upper tail event ??? Think below is correct as eventOccured has changed
for ii = 1 : nRows
    for jj = 1 : nBins        
        % Select appropriate range accounting for start  period to include the
        % initial values.
        if jj == 1
            isInRange = probForecast(ii) >= binUnitLocation( jj ) & ...
                probForecast(ii) <= binUnitLocation( jj + 1 );
            
        else
            isInRange = probForecast(ii) > binUnitLocation( jj ) & ...
                probForecast(ii) <= binUnitLocation( jj + 1 );
            
        end
                
        if isInRange
            Bin(jj).bin(end + 1, 1)  = probForecast(ii);
            Occ(jj).occ(end + 1, 1)  = eventOccured(ii);
        end
    end
end

occ = sum(eventOccured, 1) / size(eventOccured, 1);

% Initialise occk
Occk(1 : nBins)     = struct('occk', 0);
Pk(1 : nBins)       = struct('pk', 0);

for ii = 1 : nBins
    
    if size(Occ(ii).occ, 1) == 1
        Occk(ii).occk   = 0;
        
    else
        Occk(ii).occk   = sum(Occ(ii).occ, 1) / (size(Occ(ii).occ, 1) - 1);
    end
    
    if size( Bin(ii).bin, 1) == 1
        Pk(ii).pk       = 0;
        
    else
        Pk(ii).pk       = sum( Bin(ii).bin, 1) / (size( Bin(ii).bin, 1) - 1);
    end
    
end

% Initialise structure for results
Results(1 : nBins)     = struct('resolution', 0, 'reliability', 0);

for ii = 1 : nBins    
    Results(ii).resolution      = (size(Occ(ii).occ, 1) - 1)  * (Occk(ii).occk  - occ).^2;
    Results(ii).reliability     = (size(Occ(ii).occ, 1) - 1)  * (Pk(ii).pk  - Occk(ii).occk)^2;    
end

% Calculate components of the Brier score.
uncert      = occ * (1 - occ);
Res         = (1 / size(probForecast, 1)) * sum( [Results(:).resolution] );
Rel         = (1 / size(probForecast, 1)) * sum( [Results(:).reliability] );

bscore_test = uncert + Rel - Res;

BrierScore.uncertainty              = uncert;
BrierScore.resolution               = Res;
BrierScore.reliability              = Rel;
BrierScore.brierScoreDecomposition  = bscore_test;



end

