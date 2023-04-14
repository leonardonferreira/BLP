function PITs_Plot(PITs,noBin,NameModel)

mm = size(PITs,1); 


        SDbin=(1.96*sqrt(((1/noBin)*(1-1/noBin))/mm))*noBin; 
        hs = histc(PITs,0:1/noBin:1);
        bar(0:1/noBin:1,(hs./mm).*noBin,'histc')
        hold on; plot(0:1/noBin:1,ones(noBin+1,1), '-r','LineWidth',1)
        hold on; plot(0:1/noBin:1,(1+SDbin)*ones(noBin+1,1), '--r','LineWidth',1)
        hold on; plot(0:1/noBin:1,(1-SDbin)*ones(noBin+1,1), '--r','LineWidth',1)
        xlabel('Probabilities'); 
        ylabel('Density'); 
        xlim([0 1]); 
        title(NameModel);
