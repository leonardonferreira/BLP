
%hyperprior for \lambda as a function of the horizon

%-VAR---------------------------------------------------------------------%
mode=.2;    sd=.4;

x=0:.001:10;
shape=(2+mode^2/sd^2+sqrt((4+mode^2/sd^2)*mode^2/sd^2))/2;
scale=sqrt(sd^2/shape);
    
yVAR=x.^(shape-1).*exp(-x./scale)*scale^-shape/gamma(shape); 

%-variance increasing in h------------------------------------------------%
mode=.2;    stdev=.5:-.05:0.05;

x=0:.001:10;    y=nan(length(stdev),length(x));

j=1;
for sd=stdev

    shape=(2+mode^2/sd^2+sqrt((4+mode^2/sd^2)*mode^2/sd^2))/2;
    scale=sqrt(sd^2/shape);
    
    y(j,:)=x.^(shape-1).*exp(-x./scale)*scale^-shape/gamma(shape); j=j+1;

end

figure;
plot(y');
hold on
pl=plot(yVAR,'--k','LineWidth',1.5);
hold off; grid on; xlim([0 1000])
ph=legend(pl,'VAR'); set(ph,'FontSize',12)

title('Gamma(\lambda):: mode = .2, variance = .5 : 2','FontSize',12)
annotation('arrow',[.3 .3],[.6 .2]); text(150,.1,'increasing h')


% set(gcf,'PaperUnits','centimeters','PaperSize',[15 10]) %[x y]
% set(gcf,'PaperPosition',[-1 0 17 10]) %[left bottom width height]
% print(gcf,'-dpdf','GammaLambdaHfunction1.pdf'); 


%-mode increasing in h----------------------------------------------------%
modeset=.25:.05:3;
sd=.4;

x=0:.001:10;
y=nan(length(stdev),length(x));

j=1;
for mode=modeset

    shape=(2+mode^2/sd^2+sqrt((4+mode^2/sd^2)*mode^2/sd^2))/2;
    scale=sqrt(sd^2/shape);
    
    y(j,:)=x.^(shape-1).*exp(-x./scale)*scale^-shape/gamma(shape); j=j+1;

end

figure;
plot(y');
hold on
pl=plot(yVAR,'--k','LineWidth',1.5);
hold off; grid on; xlim([0 1000]); ylim([0 1.7])
ph=legend(pl,'VAR'); set(ph,'FontSize',12)

title('Gamma(\lambda):: mode = .2 : 3, variance = .4','FontSize',12)
annotation('arrow',[.4 .8],[.7 .7]); text(500,1.3,'increasing h')

% set(gcf,'PaperUnits','centimeters','PaperSize',[15 10]) %[x y]
% set(gcf,'PaperPosition',[-1 0 17 10]) %[left bottom width height]
% print(gcf,'-dpdf','GammaLambdaHfunction2.pdf'); 


