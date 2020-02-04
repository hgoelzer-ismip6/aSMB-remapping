% Plotting remapped aSMB

clear

% OBS
iism = 0;
amod = 'obs';

% MAR39
d0=load('../Data/MAR/DSMB_MARv3.9_MIROC5_rcp85.mat');
modscen='M39_MIROC5_rcp85';
% SMB mask
mask0=double(d0.MSK==5);

% Ice sheet mask
mask=mask0;
mask(isnan(mask))=0;

% Basins
load ../Data/Basins/ExtBasinScale25_nn7_50_05000m.mat wbas

% Remapping results
load ../Models/OBS/aSMB.mat

% colormap
load cmap_dsmb.mat

% load label positions
load ../Data/Basins/bc_05000m

% observed mask
nanmask0 = mask0;
nanmask0(find(nanmask0==0))=NaN;

% modelled mask
nanmask = mask;
nanmask(find(nanmask==0))=NaN;

% combined mask
nanmask12 = double(mask+mask0 > 0);
nanmask12(find(nanmask12==0))=NaN;

% 17 levels = 16 intervals is best
levels = [-6:0.4:0.4];

% Original DSMB on observed mask
nanshade_nonlinear(dsd.*nanmask0,levels)
colormap(cmap)
set(gca,'Xtick', [])
set(gca,'Ytick', [])
text(301,-18,'[m yr^{-1}]','FontSize',14,'Interpreter','tex');
hold on
contour(wbas.basinIDs',[1:25],'k')
text(bc(1,:),bc(2,:),num2str([1:25]'),'Color',[0,0,0],'Fontsize',12,'FontWeight','bold')
if(iism ~= 0)
    contour(mask0',[0.5,0.5],'Color',[0.5,0.5,0.5],'LineWidth',1.5);
end
set(gca,'color',[0.9 0.9 0.9])
axis equal
print('-dpng', '-r300', ['../Plotting/Plots/dsmb_div_' modscen '_orgobs']) 
%export_fig(['../Plotting/Plots/dsmb_div_' modscen '_orgobs'], '-png', '-r300', '-nocrop', '-transparent');

if(iism ~= 0)
% Extended DSMB on modelled mask
nanshade_nonlinear(dsd.*nanmask,levels)
colormap(cmap)
set(gca,'Xtick', [])
set(gca,'Ytick', [])
text(301,-18,'[m yr^{-1}]','FontSize',14,'Interpreter','tex');
hold on
contour(wbas.basinIDs',[1:25],'k')
text(bc(1,:),bc(2,:),num2str([1:25]'),'Color',[0,0,0],'Fontsize',12,'FontWeight','bold')
if(iism ~= 0)
    contour(mask0',[0.5,0.5],'Color',[0.5,0.5,0.5],'LineWidth',1.5)
end
set(gca,'color',[0.9 0.9 0.9])

export_fig(['../Plotting/Plots/dsmb_div_' modscen '_orgmod_' amod], '-png', '-r300', '-nocrop', '-transparent') 
end

% Remapped DSMB on modelled mask
nanshade_nonlinear(dsd_re.*nanmask,levels)
colormap(cmap)
set(gca,'Xtick', [])
set(gca,'Ytick', [])
text(301,-18,'[m yr^{-1}]','FontSize',14,'Interpreter','tex');
hold on
contour(wbas.basinIDs',[1:25],'k')
text(bc(1,:),bc(2,:),num2str([1:25]'),'Color',[0,0,0],'Fontsize',12,'FontWeight','bold')
if(iism ~= 0)
    contour(mask0',[0.5,0.5],'Color',[0.5,0.5,0.5],'LineWidth',1.5)
end
set(gca,'color',[0.9 0.9 0.9])
axis equal

print('-dpng', '-r300', ['../Plotting/Plots/dsmb_div_' modscen '_mapmod_' amod ]) 
%export_fig(['../Plotting/Plots/dsmb_div_' modscen '_mapmod_' amod ], '-png', '-r300', '-nocrop', '-transparent');

if(iism ~= 0)
%% Per basin integrals in comparison
figure
set(groot,'DefaultAxesColorOrder', parula(3))
bar(1e-9*[bint;bint_ex;bint_re]','grouped')
ax=axis;
axis([0,26, ax(3),0])
legend({'original','extended','remapped'},'Location','southeast')
ylabel('Integrated aSMB [km^3 yr^{-1}]','Interpreter','Tex')
xlabel('Basin Id')
print('-dpng', '-r300', ['../Plotting/Plots/dsmb_basinint_' modscen 'mod_map_' amod]) 
%export_fig(['../Plotting/Plots/dsmb_basinint_' modscen 'mod_map_' amod], '-png', '-r300', '-nocrop', '-transparent');
end

if(iism ~= 0)
% Mask differences 
nanshade((dsd.*mask-dsd.*mask0).*nanmask12)
load cmap_polar.mat
colormap(cmap)
caxis([-4 4])
set(gca,'Xtick', [])
set(gca,'Ytick', [])
text(301,-18,'[m yr^{-1}]','FontSize',14,'Interpreter','tex');
hold on
contour(wbas.basinIDs',[1:25],'k')
if(iism ~= 0)
    contour(mask0',[0.5,0.5],'Color',[0.5,0.5,0.5],'LineWidth',1.5);
end
text(bc(1,:),bc(2,:),num2str([1:25]'),'Color',[0,0,0],'Fontsize',12,'FontWeight','bold')
%contour(mask0',[0.5,0.5],'g')
set(gca,'color',[0.9 0.9 0.9])
axis equal
print('-dpng', '-r300', ['../Plotting/Plots/ddsmb_div_' modscen '_ext-org_' amod]) 
%export_fig(['../Plotting/Plots/ddsmb_div_' modscen '_ext-org_' amod], '-png', '-r300', '-nocrop', '-transparent');
end

% Differences of mapped to the extended DSMB
nanshade((dsd_re.*mask-dsd.*mask).*nanmask12)
load cmap_polar.mat
colormap(cmap)

caxis([-4 4])
if(iism == 0)
caxis([-2 2])
end
set(gca,'Xtick', [])
set(gca,'Ytick', [])
text(301,-18,'[m yr^{-1}]','FontSize',14,'Interpreter','tex');
hold on
contour(wbas.basinIDs',[1:25],'k')
text(bc(1,:),bc(2,:),num2str([1:25]'),'Color',[0,0,0],'Fontsize',12,'FontWeight','bold')
if(iism ~= 0)
    contour(mask0',[0.5,0.5],'Color',[0.5,0.5,0.5],'LineWidth',1.5)
end
%contour(mask0',[0.5,0.5],'g')
set(gca,'color',[0.9 0.9 0.9])
axis equal
print('-dpng', '-r300', ['../Plotting/Plots/ddsmb_div_' modscen '_map-ext_' amod]) 
%export_fig(['../Plotting/Plots/ddsmb_div_' modscen '_map-ext_' amod], '-png', '-r300', '-nocrop', '-transparent');

if(iism ~= 0)
%% Per basin differences
figure
set(groot,'DefaultAxesColorOrder', parula(2))
%bar(1e-9*[bint_ex-bint;bint_re-bint_ex]','grouped')
bar(1e-9*[bint_ex-bint;bint_ex-bint_re]','grouped')
ax=axis;
axis([0,26, ax(3),ax(4)])
%legend({'extended-original','remapped-extended'},'Location','southeast')
legend({'extended-original','extended-remapped'},'Location','southeast')
ylabel('Integrated aSMB difference [km^3 yr^{-1}]','Interpreter','Tex')
xlabel('Basin Id')
print('-dpng', '-r300', ['../Plotting/Plots/ddsmb_basinint_' modscen 'mod_map_' amod]) 
%export_fig(['../Plotting/Plots/ddsmb_basinint_' modscen 'mod_map_' amod], '-png', '-r300', '-nocrop', '-transparent');
end

if(iism == 0)
figure
set(groot,'DefaultAxesColorOrder', parula(2))
%bar(1e-9*[lookup.bint;bint_re]','grouped')
bar(1e-9*[bint;bint_re]','grouped')
tot1 = 1e-9*sum(bint);
tot2 = 1e-9*sum(bint_re);
%text(10,-70,['tot:', num2str(round(tot1))],'FontSize',14)
%text(10,-75,['tot:', num2str(round(tot2))],'FontSize',14)
axis([0,26,-80,0])
ylabel('Integrated aSMB [km^3 yr^{-1}]','Interpreter','Tex')
legend({['original (tot=', num2str(round(tot1)), ')'],['reconstructed (tot=', num2str(round(tot2)) ,')']},'Location','southwest')
xlabel('Basin Id')
print('-dpng', '-r300', ['../Plotting/Plots/dsmb_basinint_' modscen '_re']) 
%export_fig(['../Plotting/Plots/dsmb_basinint_' modscen '_re'], '-png', '-r300', '-nocrop', '-transparent');
end

if(iism ~= 0)
% elevation difference
dsur = sur-d0.SH; 
dsur(find(isnan(dsur)))=0.0;
nanshade(dsur.*nanmask12)
load cmap_polar.mat
colormap(cmap)
hold on
contour(wbas.basinIDs',[1:25],'k')
text(bc(1,:),bc(2,:),num2str([1:25]'),'Color',[0,0,0],'Fontsize',12,'FontWeight','bold')

caxis([-2000 2000])
set(gca,'Xtick', [])
set(gca,'Ytick', [])
text(345,-18,'[m]','FontSize',14,'Interpreter','tex');
set(gca,'color',[0.9 0.9 0.9])
axis equal
if(iism ~= 0)
    contour(mask0',[0.5,0.5],'Color',[0.5,0.5,0.5],'LineWidth',1.5)
end
print('-dpng', '-r300', ['../Plotting/Plots/dsur_div' modscen '_re']) 
%export_fig(['../Plotting/Plots/dsur_div' modscen '_re'], '-png', '-r300', '-nocrop', '-transparent');
end

