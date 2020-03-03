% Plot final thickness changes

clear

% Settings
modscen='M39_MIROC5-rcp85';
amod = 'OBS';

% models
t0=load(['../Data/output/transient0_' modscen '_' amod '.mat']);
t2=load(['../Data/output/transient2_' modscen '_' amod '.mat']);
t3=load(['../Data/output/transient3_' modscen '_' amod '.mat']);
t5=load(['../Data/output/transient5_' modscen '_' amod '.mat']);


load cmap_polar.mat
cmap_polar = cmap;
load cmap_dsmb.mat
cmap_dsmb = cmap;

% surface change due to DSMB
dh_c = t0.tsur(:,:,end)-t0.tsur(:,:,1);
dh_c(isnan(dh_c))=0;
dh_0 = t0.tsur_re(:,:,end)-t0.tsur_re(:,:,1);
dh_0(isnan(dh_0))=0;

% surface change due to DSMB and height feedback
dh_2 = t2.tsur_re(:,:,end)-t2.tsur_re(:,:,1);
dh_2(isnan(dh_2))=0;

% surface change due to DSMB and dsmb/dz height feedback
dh_3 = t3.tsur_re(:,:,end)-t3.tsur_re(:,:,1);
dh_3(isnan(dh_3))=0;

% surface change due to DSMB and dsmb/dz height feedback
dh_5 = t5.tsur_re(:,:,end)-t5.tsur_re(:,:,1);
dh_5(isnan(dh_5))=0;

% error in parameterisation
dh_par = t0.tsur_re(:,:,end)-t0.tsur(:,:,end);
dh_par(isnan(dh_par))=0;

% surface change
shade_nt(dh_c)
hold on
contour(dh_c',[0.1,0.1], 'Color', [0.5,0.5,0.5])
%title('sur 0')
caxis([-150 12.5])
%caxis([-150 12.5]*2.5)
colormap(cmap_dsmb)
text(350,-18,'[m]','FontSize',14,'Interpreter','tex');
print('-r300', '-dpng', ['Plots/Dh_c_2100_' modscen '_' amod]);

shade_nt(dh_0) 
hold on
contour(dh_0',[0.1,0.1], 'Color', [0.5,0.5,0.5])
%title('sur re 0')
caxis([-150 12.5])
%caxis([-150 12.5]*2.5)
colormap(cmap_dsmb)
text(350,-18,'[m]','FontSize',14,'Interpreter','tex');
print('-r300', '-dpng', ['Plots/Dh_0_2100_' modscen '_' amod]);
%shade_nt(dh_1) 
%caxis([-150 12.5])
shade_nt(dh_2) 
hold on
contour(dh_2',[0.1,0.1], 'Color', [0.5,0.5,0.5])
%title('sur re 2')
caxis([-150 12.5])
%caxis([-150 12.5]*2.5)
colormap(cmap_dsmb)
text(350,-18,'[m]','FontSize',14,'Interpreter','tex');
print('-r300', '-dpng', ['Plots/Dh_2_2100_' modscen '_' amod]);

% error in parameterisation
shade_nt(dh_par)
%title('par err: re-sur 0')
caxis([-50 50])
%caxis([-50 50]*2)
colormap(cmap_polar)
text(350,-18,'[m]','FontSize',14,'Interpreter','tex');
print('-r300', '-dpng', ['Plots/Dh_parerr_2100_' modscen '_' amod]);

shade_nt(dh_2-dh_0) 
%hold on
%contour((dh_2-dh_0)',[0.01,0.01], 'Color', [0.5,0.5,0.5])
%title('Dlapse: re2-re0')
%caxis([-5 5])
%%caxis([-5 5]*4)
%colormap(cmap_polar)
colormap(cmap_dsmb)
caxis([-150 12.5]*0.1)
text(350,-18,'[m]','FontSize',14,'Interpreter','tex');
print('-r300', '-dpng', ['Plots/Dh_lapseD_2100_' modscen '_' amod]);

shade_nt(dh_3-dh_0) 
%hold on
%contour((dh_3-dh_0)',[0.01,0.01], 'Color', [0.5,0.5,0.5])
%title('lapse: re3-re0')
caxis([-150 12.5]*0.1)
colormap(cmap_dsmb)
text(350,-18,'[m]','FontSize',14,'Interpreter','tex');
print('-r300', '-dpng', ['Plots/Dh_lapse_2100_' modscen '_' amod]);

shade_nt(dh_5-dh_0) 
%hold on
%contour((dh_5-dh_0)',[0.01,0.01], 'Color', [0.5,0.5,0.5])
%title('lapse D+: re5-re0')
caxis([-150 12.5]*0.1)
colormap(cmap_dsmb)
text(350,-18,'[m]','FontSize',14,'Interpreter','tex');
print('-r300', '-dpng', ['Plots/Dh_lapseplusD-2100_' modscen '_' amod]);

shade_nt(dh_3+dh_2-2*dh_0) 
%hold on
%contour((dh_3+dh_2-2*dh_0)',[0.01,0.01], 'Color', [0.5,0.5,0.5])
%title('lapse D+: constructed')
caxis([-150 12.5]*0.1)
colormap(cmap_dsmb)
text(350,-18,'[m]','FontSize',14,'Interpreter','tex');
print('-r300', '-dpng', ['Plots/Dh_lapseplusD_constructed_2100_' modscen '_' amod]);

