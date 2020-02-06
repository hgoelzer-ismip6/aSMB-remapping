% save basin info
clear

bas=ncload('../Data/Basins/ISMIP6Masks25_05000m.nc');

% masks
obs=ncload('../Data/RCM/zmask_05000m.nc');
mask = obs.zmask;
mask(isnan(mask)) = 0.;
omask = 1-mask;

bas.ids = 1:25;
bas.basinIDs=zeros(size(bas.basin1));
check=zeros(size(bas.basin1));
for i=1:25
    eval(['bas.basin' num2str(i) '(isnan(bas.basin' num2str(i) '))=0.;']);
    eval(['bas.basinIDs=bas.basinIDs+bas.basin' num2str(i) '*i;']);
    eval(['check = check + bas.basin' num2str(i) ';']);
end
%shade(check)

x1=1:size(bas.basinIDs,1);
y1=1:size(bas.basinIDs,2);
[y,x] = meshgrid(y1,x1);

for i=1:25
    %% mean basin position for label
    eval(['xcb=x.*(bas.basin' num2str(i) './bas.basin' num2str(i) ').*omask;']);
    eval(['ycb=y.*(bas.basin' num2str(i) './bas.basin' num2str(i) ').*omask;']);
    xcb(xcb==0) = NaN;
    ycb(ycb==0) = NaN;
    bc(2,i)=nanmedian(ycb(:));
    bc(1,i)=nanmedian(xcb(:));
end

% hand mods
bc(1,3) = bc(1,3)+17;
bc(2,3) = bc(2,3)+15;
bc(1,7) = bc(1,7)+10;
bc(2,24) = bc(2,24)+3;
bc(1,25) = bc(1,25)-15;

% make a zero basin id
bas.basinIDs(1,1)=0;
bas.basin0=bas.basin1*0;
bas.basin0(1,1) = 1;

% Plotting
shade(bas.basinIDs)
set(gca,'Xtick', [])
set(gca,'Ytick', [])
%colormap(lines(25))
hold on
contour(mask',[0.5,0.5],'k')
caxis([1,25])
text(bc(1,:),bc(2,:),num2str([1:25]'),'Color',[1,1,1],'Fontsize',12)
%print -dpng -r300 extbasins

save ../Data/Basins/ExtBasinMasks25_05000m bas

save ../Data/Basins/bc_05000m bc

