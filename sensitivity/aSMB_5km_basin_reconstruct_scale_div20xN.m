% Reconstruct aSMB based on lookup table

%clear

%mul = 5;
%iscen = 5;

% Parameters
flg_weigh = 1;
flg_match = 0;
flg_plot = 0;

colors=get(0,'DefaultAxesColorOrder');

% basin definition
load(['Data/Basins/ExtBasinMasks20x' num2str(mul) '.mat']);
x1=1:size(bas.basinIDs,1);
y1=1:size(bas.basinIDs,2);
[y,x] = meshgrid(y1,x1);
nb = 20*mul;

% area factors
da = ncload('../Data/Grid/af2_ISMIP6_GrIS_05000m.nc');
af2 = da.af2(:,:);

% dim
dx=5000;dy=5000;

% masks
obs=ncload('Data/obs1_05.nc');

% basin weights
load(['Data/Basins/ExtBasinScale20x' num2str(mul) '_div6_50'], 'wbas');

%% the original forcing file
if (iscen == 1)
    d0=load('Data/DSMB_MARv3.9_MIROC5-rcp85.mat');
    load(['lookup_MAR39_MIROC5_rcp85_b20x' num2str(mul)])
    modscen='M39_MIROC5_rcp85';
end
if (iscen == 2)
    d0=load('Data/DSMB_MARv3.9_NorESM1-rcp85.mat');
    load(['lookup_MAR39_NorESM1-rcp85_b20x' num2str(mul)])
    modscen='M39_NorESM1-rcp85';
end
if (iscen == 3)
    d0=load('Data/DSMB_MARv3.9_HadGEM2-ES-rcp85.mat');
    load(['lookup_MAR39_HadGEM2-ES-rcp85_b20x' num2str(mul)])
    modscen='M39_HadGEM2-ES-rcp85';
end
if (iscen == 4)
    d0=load('Data/DSMB_MARv3.9_IPSL-CM5-MR-rcp85.mat');
    load(['lookup_MAR39_IPSL-CM5-MR-rcp85_b20x' num2str(mul)])
    modscen='M39_IPSL-CM5-MR-rcp85';
end
if (iscen == 5)
    d0=load('Data/DSMB_MARv3.9_CSIRO-Mk3.6-rcp85.mat');
    load(['lookup_MAR39_CSIRO-Mk36-rcp85_b20x' num2str(mul)])
    modscen='M39_CSIRO-Mk36-rcp85';
end
if (iscen == 6)
    d0=load('Data/DSMB_MARv3.9_ACCESS1.3-rcp85.mat');
    load(['lookup_MAR39_ACCESS13-rcp85_b20x' num2str(mul)])
    modscen='M39_ACCESS13-rcp85';
end


% SMB product mask
mask0=double(d0.MSK==5);

%% Load a modelled geometry for reconstruction
amod = 'obs';
nc = obs;
nc2 = obs;
nc2.sftgif= mask0;
sur=nc.orog;

% Ice sheet mask
mask=nc2.sftgif;
mask(isnan(mask))=0;

lat=d0.LAT;
dsd=d0.DSMB;
dsd_re=zeros(size(dsd));

bint_re=zeros(1,nb);
bint_ex=zeros(1,nb);
bint=zeros(1,nb);

% loop through basins
for i=1:nb

    %% set current basin, mask and lookup
    eval(['sur_b=sur.*(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);
    eval(['mask_b=(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);
    eval(['ima_b=mask.*(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);

    eval(['look0=lookup.b' num2str(wbas.n0(i)) ';']);
    %% set neighbor basin and lookup
    eval(['look1=lookup.b' num2str(wbas.n1(i)) ';']);
    eval(['look2=lookup.b' num2str(wbas.n2(i)) ';']);
    eval(['look3=lookup.b' num2str(wbas.n3(i)) ';']);
    eval(['look4=lookup.b' num2str(wbas.n4(i)) ';']);
    eval(['look5=lookup.b' num2str(wbas.n5(i)) ';']);
    eval(['look6=lookup.b' num2str(wbas.n6(i)) ';']);
    
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
rmsn = nanrms(dd(:));
[iscen, mul, rmsn ]

if(flg_plot)
% Plotting
load cmap_dsmb.mat

% load label positions
load(['bc20x' num2str(mul) '.mat'],'bc')

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
contour(wbas.basinIDs',[1:nb],'k')
text(bc(1,:),bc(2,:),num2str([1:nb]'),'Color',[0,0,0],'Fontsize',12)
contour(mask0',[0.5,0.5],'r')

print('-dpng', '-r300', ['../Plotting/Plots/dsmb_div_' modscen '_orgobs' num2str(mul)]) 

% Extended DSMB on modelled mask
nanshade_nonlinear(dsd.*nanmask,levels)
colormap(cmap)
set(gca,'Xtick', [])
set(gca,'Ytick', [])
text(301,-18,'[m yr^{-1}]','FontSize',14,'Interpreter','tex');
hold on
contour(wbas.basinIDs',[1:nb],'k')
text(bc(1,:),bc(2,:),num2str([1:nb]'),'Color',[0,0,0],'Fontsize',12)
contour(mask0',[0.5,0.5],'r')

print('-dpng', '-r300', ['../Plotting/Plots/dsmb_div_' modscen '_orgmod_' amod num2str(mul)]) 


% Remapped DSMB on modelled mask
nanshade_nonlinear(dsd_re.*nanmask,levels)
colormap(cmap)
set(gca,'Xtick', [])
set(gca,'Ytick', [])
text(301,-18,'[m yr^{-1}]','FontSize',14,'Interpreter','tex');
hold on
contour(wbas.basinIDs',[1:nb],'k')
text(bc(1,:),bc(2,:),num2str([1:nb]'),'Color',[0,0,0],'Fontsize',12)
contour(mask0',[0.5,0.5],'r')

print('-dpng', '-r300', ['../Plotting/Plots/dsmb_div_' modscen '_mapmod_' amod '_wgt' num2str(flg_weigh) num2str(mul)]) 


%% Per basin integrals in comparison
figure
bar(1e-9*[bint;bint_ex;bint_re]','grouped')
ax=axis;
axis([0,nb+1, ax(3),0])
legend({'original','extended','remapped'},'Location','southeast')
ylabel('Integrated aSMB [km^3]')
xlabel('Basin Id')
print('-dpng', '-r300', ['../Plotting/Plots/dsmb_basinint_' modscen 'mod_map_' amod '_wgt' num2str(flg_weigh) num2str(mul)]) 


% Mask differences 
nanshade((dsd.*mask-dsd.*mask0).*nanmask12)
load cmap_polar.mat
colormap(cmap)
caxis([-4 4])
set(gca,'Xtick', [])
set(gca,'Ytick', [])
text(301,-18,'[m yr^{-1}]','FontSize',14,'Interpreter','tex');
hold on
contour(wbas.basinIDs',[1:nb],'k')
text(bc(1,:),bc(2,:),num2str([1:nb]'),'Color',[0,0,0],'Fontsize',12)
%contour(mask0',[0.5,0.5],'g')
print('-dpng', '-r300', ['../Plotting/Plots/ddsmb_div_' modscen '_ext-org_' amod '_wgt' num2str(flg_weigh) num2str(mul)]) 

% Differences of mapped to the extended DSMB
nanshade((dsd_re.*mask-dsd.*mask).*nanmask12)
load cmap_polar.mat
colormap(cmap)
caxis([-4 4])
set(gca,'Xtick', [])
set(gca,'Ytick', [])
text(301,-18,'[m yr^{-1}]','FontSize',14,'Interpreter','tex');
hold on
contour(wbas.basinIDs',[1:nb],'k')
text(bc(1,:),bc(2,:),num2str([1:nb]'),'Color',[0,0,0],'Fontsize',12)
%contour(mask0',[0.5,0.5],'g')
print('-dpng', '-r300', ['../Plotting/Plots/ddsmb_div_' modscen '_map-ext_' amod '_wgt' num2str(flg_weigh) num2str(mul)]) 

%% Per basin differences
figure
bar(1e-9*[bint_ex-bint;bint_re-bint_ex]','grouped')
ax=axis;
axis([0,nb+1, ax(3),ax(4)])
legend({'extended-original','remapped-extended'},'Location','southeast')
ylabel('Integrated aSMB difference [km^3]','Interpreter','Tex')
xlabel('Basin Id')
print('-dpng', '-r300', ['../Plotting/Plots/ddsmb_basinint_' modscen 'mod_map_' amod '_wgt' num2str(flg_weigh) num2str(mul)]) 


figure
%bar(1e-9*[lookup.bint;bint_re]','grouped')
bar(1e-9*[bint;bint_re]','grouped')
ax=axis;
axis([0,nb+1, ax(3),ax(4)])
ylabel('Integrated DSMB [km^3]','Interpreter','Tex')
legend({'original','reconstructed'},'Location','southeast')
xlabel('Basin Id')
print('-dpng', '-r300', ['../Plotting/Plots/dsmb_basinint_' modscen '_re' num2str(mul)]) 

end
