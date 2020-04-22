% Plot basins and weighting
clear

% load obs
dg = ncload('../Data/RCM/grid_MARv3.9_05000m.nc');
dm = ncload('../Data/Grid/zmask_05000m.nc');

cma = (dg.MSK<2);

load ../Data/Basins/ExtBasinScale25_nn7_50_05000m wbas

% masking
mask0=double(dg.MSK==5);


% label positions
load ../Data/Basins/bc_05000m.mat


% basin weight
%shade(wbas.wg)
shade_nt(wbas.wgc0)
hold on
contour(wbas.basinIDs',[1:1:25],'k','Linewidth',1)
caxis([0.25 1])
contour(cma',[0.5,0.5],'Color','r','Linewidth',2)
contour(mask0',[0.5,0.5],'Linewidth',1,'Color','w','Linewidth',2)
text(bc(1,:),bc(2,:),num2str([1:25]'),'Color',[0,0,0],'Fontsize',12,'FontWeight','bold')
set(gca,'XTicklabels',[])
set(gca,'YTicklabels',[])
set(gca,'XTick',[])
set(gca,'YTick',[])

%print -dpng -r300 Plots/basinWeights_div6_50
export_fig(['Plots/basinWeights_div6_50'], '-png', '-r600', '-nocrop', '-transparent');

c = [0,0,0];
set(groot,'DefaultAxesColorOrder', c)
shade_nb(mask0)
c = colormap(gray);
colormap(c)
hold on
contour(wbas.basinIDs',[0:25],'Color',[0.5,0.5,0.5],'Linewidth',1)
contour(cma',[0.5,0.5],'Color',[0.2,0.2,0.2],'Linewidth',1)
text(bc(1,:),bc(2,:),num2str([1:25]'),'Color',[0,0,0],'Fontsize',13,'FontWeight','bold')
set(gca,'XTicklabels',[])
set(gca,'YTicklabels',[])
set(gca,'XTick',[])
set(gca,'YTick',[])

%caxis([-3,1])
caxis([-5,1])

print -dpng -r300 Plots/basinIDs_div
export_fig(['Plots/basinIDs_div'], '-png', '-r600', '-nocrop', '-transparent');
