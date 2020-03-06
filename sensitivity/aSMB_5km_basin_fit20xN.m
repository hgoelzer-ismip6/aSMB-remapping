% Construct lookup tables for remapping

%clear

%mul = 5;
%iscen = 4;

% fill with nan for high elevation?
flg_nanfill = 1;
flg_plot = 0;

colors=get(0,'DefaultAxesColorOrder');

% basin definition
load(['Data/Basins/ExtBasinMasks20x' num2str(mul) '.mat']);

% area factors
da = ncload('../Data/Grid/af2_ISMIP6_GrIS_05000m.nc');
af2 = da.af2(:,:);

% dim
dx=5000;dy=5000;

% number of basins
nb = length(bas.ids);

% masks
obs=ncload('Data/obs1_05.nc');


% ISMIP6 MAR3.9 - conservative interpolated ref1960-1998
if (iscen == 1)
d0=load('Data/DSMB_MARv3.9_MIROC5-rcp85.mat');
lookup_file=['lookup_MAR39_MIROC5_rcp85_b20x' num2str(mul)];
sur=d0.SH;
end
if (iscen == 2)
d0=load('Data/DSMB_MARv3.9_NorESM1-rcp85.mat');
lookup_file=['lookup_MAR39_NorESM1-rcp85_b20x' num2str(mul)];
sur=d0.SH;
end
if (iscen == 3)
d0=load('Data/DSMB_MARv3.9_HadGEM2-ES-rcp85.mat');
lookup_file=['lookup_MAR39_HadGEM2-ES-rcp85_b20x' num2str(mul)];
sur=d0.SH;
end
if (iscen == 4)
d0=load('Data/DSMB_MARv3.9_IPSL-CM5-MR-rcp85.mat');
lookup_file=['lookup_MAR39_IPSL-CM5-MR-rcp85_b20x' num2str(mul)];
sur=d0.SH;
end
if (iscen == 5)
d0=load('Data/DSMB_MARv3.9_CSIRO-Mk3.6-rcp85.mat');
lookup_file=['lookup_MAR39_CSIRO-Mk36-rcp85_b20x' num2str(mul)];
sur=d0.SH;
end
if (iscen == 6)
d0=load('Data/DSMB_MARv3.9_ACCESS1.3-rcp85.mat');
lookup_file=['lookup_MAR39_ACCESS13-rcp85_b20x' num2str(mul)];
sur=d0.SH;
end


%mask=d0.MSK > 0;
mask=obs.zmask;
dsd=d0.DSMB.*(mask./double(mask));
%dsd=d0.DSMB.*(mask);
bint=zeros(1,nb);

f = figure;
set(f, 'DefaultLineLineWidth', 0.5)
set(f, 'DefaultLineMarkerSize', 4)
set(0, 'DefaultAxesFontSize', 8)

% elevation levels
sstep = 100;
ds = 100;
ss=ds:sstep:3500;

[p,n] = numSubplots(nb);

oldlook = zeros(2,length(ss)+1); 

for i=1:nb
%    figure
    subplot(p(2),p(1),i)
    hold on; box on;
    eval(['dsd_b=dsd.*(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);
    eval(['sur_b=sur.*(bas.basin' num2str(i) './bas.basin' num2str(i) ');']);
    
    plot(sur_b(:),dsd_b(:),'.');
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

        plot(ssel,ysel,'.','Color',colors(n,:))

        %% fill nans with last value
        if(isnan(y1) && flg_nanfill)
            look(:,n)=[s0,yold];
            plot(s0,yold,'o','Color',colors(n,:))
        else
            look(:,n)=[s0,y1];
            yold=y1;
            %plot(s0,y1,'x','Color',colors(n,:))
        end
        %% fill first value for h=0 from second
        look(2,1) = look(2,2);

    end

    %% if first values are zero, set neighbor basin
    lastmatch = 0;
    for k = 1:20 % search up 2000 m
        if look(2,k) == 0.;
            lastmatch = k;
        else
            break;
        end
    end
    look(2,1:lastmatch) = oldlook(2,1:lastmatch); 
    % remember
    oldlook = look;
    
    plot(look(1,:),look(2,:),'-k')
    
    eval(['lookup.b' num2str(i) ' = look;']);
    axis([0 3300 -5 1])
%    set(gca,'XTick',[])
%    axis([0 3300 -2 1])
%    title(['B' num2str(bas.ids(i)) ' ID' num2str(i) ])
    text(2300,-4,['ID ', num2str(bas.ids(i))])    
end

% add dummy lookup
look0 = look;
look0(2,:) = 1:length(look0);
lookup.b0 = look0;

%% add basin totals to lookup
lookup.bint = bint;

%% write lookup table for model/scenario
save(lookup_file,'lookup')

% print('-dpng', '-r300', ['../Plotting/Plots/bfit_gradients_sd20x' num2str(mul)])
