% Construct lookup tables for remapping

clear

% fill with nan for high elevation?
flg_nanfill = 1;
flg_plot = 1;

% basin definition
load ../Data/Basins/ExtBasinMasks25_05000m.mat

% area factors
da = ncload('../Data/Grid/af2_ISMIP6_GrIS_05000m.nc');
af2 = da.af2(:,:);

% dim
dx=5000;dy=5000;

% number of basins
nb = length(bas.ids);

% masks
obs=ncload('../Data/RCM/zmask_05000m.nc');
mask=obs.zmask;

% 0=initMIP; 1=MIROC8.5; 2=NorESM8.5; 3=CANSM8.5; 4=MIROC4.5; 5=M37 MIROC8.5;
% 6=MAR39 MIROC8.5
iscen = 6;

if (iscen ==0) % initmip
d0=ncload('../Data/initMIP/dsmb_05e3413_ISMIP6_v2.nc');
lookup_file='lookup_initMIP_b25';
obs=ncload('../Data/RCM/orog_05000m.nc');
sur = obs.orog;
end

if (iscen == 1)
d0=load('../Data/RCM/DSMB_MARv3.5_MIROC5_rcp85.mat');
lookup_file='lookup_MIROC5_rcp85_b25';
sur=d0.SH;
end
if (iscen == 2)
d0=load('../DATA/RCM/DSMB_MARv3.5_NorESM1_rcp85.mat');
lookup_file='lookup_NorESM1_rcp85_b25';
sur=d0.SH;
end
if (iscen == 3)
d0=load('../DATA/RCM/DSMB_MARv3.5_CanESM2_rcp85.mat');
lookup_file='lookup_CanESM2_rcp85_b25';
sur=d0.SH;
end
if (iscen == 4)
d0=load('../DATA/RCM/DSMB_MARv3.5_MIROC5_rcp45.mat');
lookup_file='lookup_MIROC5_rcp45_b25';
sur=d0.SH;
end
% new data
if (iscen == 5)
d0=load('../Data/RCM/DSMB_MARv3.7_MIROC5_rcp85.mat');
lookup_file='lookup_MAR37_MIROC5_rcp85_b25';
sur=d0.SH;
end
if (iscen == 6)
d0=load('../Data/RCM/DSMB_MARv3.9_MIROC5_rcp85.mat');
lookup_file='../Data/lookup/aSMB_lookup_b25_MARv3.9-MIROC5-rcp85.mat';
sur=d0.SH;
end

%mask=d0.MSK > 0;
dsd=d0.DSMB.*(mask./double(mask));
%dsd=d0.DSMB.*(mask);
bint=zeros(1,25);

if (flg_plot)
    f = figure;
    set(f, 'DefaultLineLineWidth', 0.5)
    set(f, 'DefaultLineMarkerSize', 4)
    set(0, 'DefaultAxesFontSize', 8)
    colors = distinguishable_colors(36);
end

% elevation levels
sstep = 100;
ds = 100;
ss=ds:sstep:3500;

for i=1:25

    eval(['dsd_b=dsd.*(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);
    eval(['sur_b=sur.*(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);

if (flg_plot)
%    figure
    subplot(5,5,i)
    hold on; box on;
    plot(sur_b(:),dsd_b(:),'.');
end
    
    %% integral dsmb for this basin
    bint(i)=nansum(nansum(dsd_b.*af2))*dx*dy;
    
    %% fit a lookup table to the data
    %% centers at sstep intervals with ds range
    look = zeros(2,length(ss)+1); 
    n=1;
    yold=0.;
    for s0=ss
        n=n+1;
        %% find local average DSMB for given elevation range
%        [y1,ysel,s1,ssel]= find_local_average(sur_b,dsd_b,s0,ds);
        [y1,ysel,s1,ssel]= find_local_median(sur_b,dsd_b,s0,ds);

if (flg_plot)
        plot(ssel,ysel,'.','Color',colors(n,:))
end
        %% fill nans with last value
        if(isnan(y1) && flg_nanfill)
            look(:,n)=[s0,yold];
if (flg_plot)
            plot(s0,yold,'o','Color',colors(n,:))
end
        else
            look(:,n)=[s0,y1];
            yold=y1;
            %plot(s0,y1,'o','Color',colors(n,:))
        end
        %% fill first value for h=0 from second
        look(2,1) = look(2,2);

    end

    eval(['lookup.b' num2str(i) ' = look;']);
    
    if (flg_plot)
        plot(look(1,:),look(2,:),'-k','LineWidth',1)
        axis([0 3300 -5 1])
        %    set(gca,'XTick',[])
        %    axis([0 3300 -2 1])
        %    title(['B' num2str(bas.ids(i)) ' ID' num2str(i) ])
        set(gca,'FontSize',20)
        text(2500,-4,['b', num2str(bas.ids(i)) ''],'FontSize',20,'Interpreter','Tex')    
    end

end

%% add basin totals to lookup
lookup.bint = bint;

%% write lookup table for model/scenario
save(lookup_file,'lookup')

%% Write netcdf
zd = [0, ss];
nz = length(zd);
table = zeros(nz,nb);
for i=1:nb
    eval(['table(:,i) = lookup.b' num2str(i) '(2,:);']);
end
%ncwrite_DSMBTable(['Public/DSMB_nb25_' lookup_file '.nc'],zd,table,'DSMB',{'z','b'});

if (flg_plot)
    xlabel('Surface elevation [m]')
    ylabel('aSMB [m yr^{-1}]','Interpreter', 'tex')
    %print('-dpng', '-r300', ['../Plotting/Plots/bfit_gradients'])
end
