function plotsofuniformity(pdens,figtitle)

bin = 10;
m = length(pdens);

% histograms
fn = figure;
hs = histc(pdens,0:1/bin:1);
bar(0:1/bin:1,bin*hs./m,'histc')
hold on
plot(0:1/bin:1,bin*(1/bin)*ones(bin+1,1), '-r','LineWidth',2)
hold on
plot(0:1/bin:1,bin*(1/bin+1.96*sqrt((1/bin*(1-1/bin))/m))*ones(bin+1,1), '--r','LineWidth',2)
hold on
plot(0:1/bin:1,bin*(1/bin-1.96*sqrt((1/bin*(1-1/bin))/m))*ones(bin+1,1), '--r','LineWidth',2)
xlim([0 1])
title(figtitle)
box on

% tests
rvec = [0:0.001:1];
for r = 1:size(rvec,2)
    ecdf(:,r) = mean(pdens < rvec(:,r));
end
fn = figure; 
xlabel('r') 
ylabel('\phi_P(r)') 
plot(rvec,ecdf,'LineWidth',2)
hold on
plot(rvec,rvec,'r','LineWidth',2);
hold on
plot(rvec,rvec + 1.34/sqrt(m),'r:','LineWidth',2);
hold on
plot(rvec,rvec - 1.34/sqrt(m),'r:','LineWidth',2);
hold off
xlim([0 1])
ylim([0 1])
grid on
legend('Empirical','Theoretical','5% Critical Value','Location','NorthWest'); 
title(figtitle)
box on