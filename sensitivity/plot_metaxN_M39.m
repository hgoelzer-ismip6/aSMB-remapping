% plot improvement with N

clear

load meta_rms_M39
rmss_M39 = rmss;

n = [20 40 60 80 100];

n_ism = 25;
rmsn_ism = 0.0938;

mscen_M39 = {'MIROC5_rcp85', 'NorESM1-rcp85', 'HadGEM2-ES-rcp85', 'IPSL-CM5-MR-rcp85','CSIRO-Mk36-rcp85','ACCESS13-rcp85'};


% Plot settings
set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFIgurePaperPositionMode','auto')
set(0, 'DefaultAxesFontSize', 16)
set(0, 'DefaultLineLineWidth', 2)

figure
set(gcf, 'DefaultLineLineWidth', 2)
hold on; box on;
h = [];

h(1)=plot(n,rmss_M39(1,:),'--dg','MarkerSize',10);
h(2)=plot(n,rmss_M39(2,:),'--ok','MarkerSize',10);
h(3)=plot(n,rmss_M39(3,:),'--xb','MarkerSize',10);
h(4)=plot(n,rmss_M39(4,:),'--+','MarkerSize',10,'Color',[0.7,0.5,0.2]);
h(5)=plot(n,rmss_M39(5,:),'--vc','MarkerSize',10);
h(6)=plot(n,rmss_M39(6,:),'--^m','MarkerSize',10);

%h(7)=plot(n_ism,rmsn_ism,'or','MarkerSize',10,'MarkerFaceColor','r');

set(gca,'Xtick', [20,25,40,60,80,100])
ylabel('RMSE (aSMB) [m yr^{-1}]','Interpreter','Tex')
xlabel('Number of Basins')
%order=[1,2,3,4,5,6];
order=[3,1,4,6,2,5];
%order=[7,5,6,3,2,1,4];

axis([10 110 0.05 0.13])
legend(h(order),mscen_M39([order]),'Interpreter','None','FontSize',10,'Location','sw')
print('-dpng','-r300','../Plotting/Plots/metaxN_noleg_M39')
