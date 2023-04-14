function [PITs] = PITs_CalculateAndPlot(TARGET,ForecastDensities,noBin,WhatPlot)

mm = size(TARGET,1); 
PITs = zeros(mm,1);
for tt=1:mm
    PITs(tt,1) = sum(ForecastDensities(tt,:)<=TARGET(tt,1))/size(ForecastDensities,2);
end

if WhatPlot==1
        SDbin=(1.96*sqrt(((1/noBin)*(1-1/noBin))/mm))*noBin; 
        hs = histc(PITs,0:1/noBin:1);
        bar(0:1/noBin:1,(hs./mm).*noBin,'histc')
        hold on; plot(0:1/noBin:1,ones(noBin+1,1), '-r','LineWidth',1)
        hold on; plot(0:1/noBin:1,(1+SDbin)*ones(noBin+1,1), '--r','LineWidth',1)
        hold on; plot(0:1/noBin:1,(1-SDbin)*ones(noBin+1,1), '--r','LineWidth',1)
        xlabel('Probabilities'); 
        ylabel('Density'); 
        xlim([0 1]); 
end


