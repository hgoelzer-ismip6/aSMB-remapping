% Reconstruct aSMB based on lookup table

clear

addpath('../toolbox')

% Plot settings
set(0,'DefaultTextInterpreter','none');
set(groot,'DefaultFIgurePaperPositionMode','auto')
set(groot,'DefaultAxesFontSize', 16)
set(groot,'DefaultLineLineWidth', 4)

% Parameters
flg_weigh = 1;
flg_match = 0;
flg_plot = 1;

% 0=obs; 1=vub; 2=MPI; 3=JPL; 4=BGC;
%iism = 0;
iism = 1;

colors=get(0,'DefaultAxesColorOrder');

% basin definition
load ../Data/Basins/ExtBasinMasks25_05000m.mat
x1=1:size(bas.basinIDs,1);
y1=1:size(bas.basinIDs,2);
[y,x] = meshgrid(y1,x1);

% area factors
da = ncload('../Data/Grid/af2_ISMIP6_GrIS_05000m.nc');
af2 = double(da.af2(:,:));

% dim
dx=5000;dy=5000;

% param
secpyear = 31556926;

% basin weights
load ../Data/Basins/ExtBasinScale25_nn7_50_05000m.mat wbas

%% load the original forcing file
% ISMIP6 forcing
dm=ncload('../Data/RCM/aSMB_MARv3.9-yearly-MIROC5-rcp85_ltm2091-2100_e05000m.nc');
load ../Data/lookup/aSMB_lookup_b25_MARv3.9-MIROC5-rcp85.mat
modscen='M39_MIROC5-rcp85';
dg = ncload(['../Data/RCM/grid_MAR3.9_05000m.nc']);
sur=dg.SRF;
d0.DSMB = dm.aSMB * secpyear / 1000;
d0.MSK = dg.MSK;
d0.LAT = dg.LAT;
d0.SH = dg.SRF;

% SMB product mask
mask0=double(d0.MSK==5);

%% Load a modelled geometry for reconstruction
if (iism == 0)
amod = 'obs';
nc = ncload('../Data/RCM/orog_05000m.nc');
nc2.sftgif= mask0;
end
if (iism == 1)
amod = 'vubgism';
nc=ncload('../Models/VUBGISM/orog_05000m.nc');
nc2=ncload('../Models/VUBGISM/sftgif_05000m.nc');
end
if (iism == 2)
amod = 'mpimpism';
nc=ncload('../Models/MPIMPISM/orog_05000m.nc');
nc2=ncload('../Models/MPIMPISM/sftgif_05000m.nc');
end
if (iism == 3)
amod = 'jplissm';
nc=ncload('../Models/JPLISSM/orog_05000m.nc');
nc2=ncload('../Models/JPLISSM/sftgif_05000m.nc');
% correct mask problem
nc2.sftgif = double(nc2.sftgif>0);
end
if (iism == 4)
amod = 'bgcbi';
nc=ncload('../Models/BGCBISICLES/orog_05000m.nc');
nc2=ncload('../Models/BGCBISICLES/sftgif_05000m.nc');
end

sur=nc.orog;

% dummy lookup for zero
dummy0 = lookup.b1;

% Ice sheet mask
mask=nc2.sftgif;
mask(isnan(mask))=0;

lat=d0.LAT;
dsd=d0.DSMB;
dsd_re=zeros(size(dsd));
bc=zeros(2,25);

bint_re=zeros(1,25);
bint_ex=zeros(1,25);
bint=zeros(1,25);

% loop through basins
for i=1:25

    %% set current basin, mask and lookup
    eval(['sur_b=sur.*(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);
    eval(['mask_b=(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);
    eval(['ima_b=mask.*(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);

    %% set neighbor basin and lookup
    look0 = dummy0;
    if (wbas.n0(i)>0)
        eval(['look0=lookup.b' num2str(wbas.n0(i)) ';']);
    end
    look1 = dummy0;
    if (wbas.n1(i)>0)
        eval(['look1=lookup.b' num2str(wbas.n1(i)) ';']);
    end
    look2 = dummy0;
    if (wbas.n2(i)>0)
        eval(['look2=lookup.b' num2str(wbas.n2(i)) ';']);
    end
    look3 = dummy0;
    if (wbas.n3(i)>0)
        eval(['look3=lookup.b' num2str(wbas.n3(i)) ';']);
    end
    look4 = dummy0;
    if (wbas.n4(i)>0)
        eval(['look4=lookup.b' num2str(wbas.n4(i)) ';']);
    end
    look5 = dummy0;
    if (wbas.n5(i)>0)
        eval(['look5=lookup.b' num2str(wbas.n5(i)) ';']);
    end
    look6 = dummy0;
    if (wbas.n6(i)>0)
        eval(['look6=lookup.b' num2str(wbas.n6(i)) ';']);
    end

    %% use lookup table to determine DSMB
    dsd_b0 = interp1(look0(1,:),look0(2,:),sur_b);
    dsd_b1 = interp1(look1(1,:),look1(2,:),sur_b);
    dsd_b2 = interp1(look2(1,:),look2(2,:),sur_b);
    dsd_b3 = interp1(look3(1,:),look3(2,:),sur_b);
    dsd_b4 = interp1(look4(1,:),look4(2,:),sur_b);
    dsd_b5 = interp1(look5(1,:),look5(2,:),sur_b);
    dsd_b6 = interp1(look6(1,:),look6(2,:),sur_b);

    if (flg_weigh == 0)
        %% combine unweighted
        dsd_b = dsd_b0.*wbas.wgn0;
    else
        %% combine according to weights
        dsd_b = dsd_b0.*wbas.wgc0 + dsd_b1.*wbas.wgc1 + dsd_b2.*wbas.wgc2 + dsd_b3.*wbas.wgc3 + dsd_b4.*wbas.wgc4 + dsd_b5.*wbas.wgc5 + dsd_b6.*wbas.wgc6;
    end
%    shade(dsd_b)

    %% extended integral dsmb for this basin
    dsd_ex = dsd.*ima_b.*mask;
    bint_ex(i)=nansum(nansum(dsd_ex.*af2))*dx*dy;
    
    %% mask
    dsd_b = dsd_b.*mask;

    %% integral dsmb for this basin
    bint(i)=nansum(nansum(dsd.*mask_b.*mask0.*af2))*dx*dy;
    bint_re(i)=nansum(nansum(dsd_b.*af2))*dx*dy;

    if (flg_match == 1)        
        %% adjust to match total
        dsd_b = dsd_b * (bint(i)/bint_re(i));
    end

    %% check integral again 
    bint_out(i)=nansum(nansum(dsd_b.*af2))*dx*dy;

    %% replace nan by zeros to add all basins together
    dsd_b(isnan(dsd_b))=0;
    dsd_re = dsd_re+dsd_b;

end

dd = (dsd_re-dsd).*mask;
rms = nanrms(dd(:));

if (iism == 0)
    save ../Models/OBS/aSMB.mat dsd dsd_re bint bint_re bint_ex
end

if(flg_plot)
% Plotting
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
%print('-dpng', '-r300', ['Plots/dsmb_div_' modscen '_orgobs']) 
export_fig(['../Plotting/Plots/dsmb_div_' modscen '_orgobs'], '-png', '-r300', '-nocrop', '-transparent');

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

%print('-dpng', '-r300', ['Plots/dsmb_div_' modscen '_mapmod_' amod '_wgt' num2str(flg_weigh)]) 
export_fig(['../Plotting/Plots/dsmb_div_' modscen '_mapmod_' amod '_wgt' num2str(flg_weigh)], '-png', '-r300', '-nocrop', '-transparent');

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
%print('-dpng', '-r300', ['Plots/dsmb_basinint_' modscen 'mod_map_' amod '_wgt' num2str(flg_weigh)]) 
export_fig(['../Plotting/Plots/dsmb_basinint_' modscen 'mod_map_' amod '_wgt' num2str(flg_weigh)], '-png', '-r300', '-nocrop', '-transparent');
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
%print('-dpng', '-r300', ['Plots/ddsmb_div_' modscen '_ext-org_' amod '_wgt' num2str(flg_weigh)]) 
export_fig(['../Plotting/Plots/ddsmb_div_' modscen '_ext-org_' amod '_wgt' num2str(flg_weigh)], '-png', '-r300', '-nocrop', '-transparent');
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
%print('-dpng', '-r300', ['Plots/ddsmb_div_' modscen '_map-ext_' amod '_wgt' num2str(flg_weigh)]) 
export_fig(['../Plotting/Plots/ddsmb_div_' modscen '_map-ext_' amod '_wgt' num2str(flg_weigh)], '-png', '-r300', '-nocrop', '-transparent');

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
%print('-dpng', '-r300', ['Plots/ddsmb_basinint_' modscen 'mod_map_' amod '_wgt' num2str(flg_weigh)]) 
export_fig(['../Plotting/Plots/ddsmb_basinint_' modscen 'mod_map_' amod '_wgt' num2str(flg_weigh)], '-png', '-r300', '-nocrop', '-transparent');
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
axis([0,26,-100,0])
ylabel('Integrated aSMB [km^3 yr^{-1}]','Interpreter','Tex')
legend({['original (tot=', num2str(round(tot1)), ')'],['reconstructed (tot=', num2str(round(tot2)) ,')']},'Location','southwest')
xlabel('Basin Id')
%print('-dpng', '-r300', ['Plots/dsmb_basinint_' modscen '_re']) 
export_fig(['../Plotting/Plots/dsmb_basinint_' modscen '_re'], '-png', '-r300', '-nocrop', '-transparent');
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
%print('-dpng', '-r300', ['Plots/dsur_div' modscen '_re']) 
export_fig(['../Plotting/Plots/dsur_div' modscen '_re'], '-png', '-r300', '-nocrop', '-transparent');
end

end % plt_flg
