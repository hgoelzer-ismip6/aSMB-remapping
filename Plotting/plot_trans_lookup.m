% Plot the lookup tables

clear

set(groot,'DefaultFIgurePaperPositionMode','auto')
set(groot,'DefaultAxesFontSize', 16)
set(groot,'DefaultLineLineWidth', 4)

% param
secpyear = 31556926;

lookup_file='../Data/lookup/TaSMB_trans_lookup_b25_MARv3.9-MIROC5-rcp85.nc';

% basin definition
load ../Data/Basins/ExtBasinMasks25_05000m.mat

figure

% produce custom line colors
cmap = colormap(jet(17));
colororder = cmap(1:1:17,:);

set(0,'DefaultAxesColorOrder', colororder);

dl = ncload(lookup_file);
lookup.table = dl.aSMB_ltbl;
lookup.ss = dl.z;

for b=1:25
    subplot(5,5,b)
    for t=1:5:85
        hold on; box on;
        eval(['look = lookup.table(:,b,t);']);
        plot(lookup.ss,look(:) * secpyear / 1000,'LineStyle','-','Linewidth',1)
%    title(['B' num2str(bas.ids(b)) ' ID' num2str(b) ])
        title(['B' num2str(bas.ids(b))])
%    axis([0 3300 -5 2])
        axis([0 3300 -6 1])
        caxis([0 100])
    end
end

xlabel('Surface elevation [m]')
ylabel('aSMB [m yr^{-1}]','Interpreter','tex')
caxis([0 85])


% print('-dpng', '-r300', ['Plots/trans_lookup_MAR39_MIROC5_rcp85_b25.png']); 

% cb = colorbar;
% set(cb,'Ticks',[5:15:85])


% print('-dpng', '-r300', ['Plots/trans_lookup_MAR39_MIROC5_rcp85_b25_leg.png']); 
