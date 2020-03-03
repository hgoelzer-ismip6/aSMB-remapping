% Run transient model adjusting elevation as we go
% call with meta_trans_dsmb.m
%clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
flg_weigh = 1;
flg_match = 0;

% flag for experiment type
%flg_t=0;
%flg_t=2;
%flg_t=3;
%flg_t=5;

% flag for plotting 
flg_plot=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load cmap_dsmb

colors=get(0,'DefaultAxesColorOrder');

% basin definition
load(['../Data/Basins/ExtBasinMasks25_05000m.mat']);
x1 = 1:size(bas.basinIDs,1);
y1 = 1:size(bas.basinIDs,2);
nb = length(bas.ids);
[y,x] = meshgrid(y1,x1);

% area factors
da = ncload('../Data/Grid/af2_ISMIP6_GrIS_05000m.nc');
af2 = double(da.af2(:,:));

% dim
dx=5000;dy=5000;

% param
secpyear = 31556926;
rhof = 1000;
rhoi = 910;

% basin weights
load(['../Data/Basins/ExtBasinScale25_nn7_50_05000m.mat'], 'wbas');

% ISMIP6 forcing
d0=ncload(['../Data/RCM/aSMB_MARv3.9-yearly-MIROC5-rcp85-2015-2100_05000m.nc']);
d1=ncload(['../Data/RCM/dSMBdz_MARv3.9-yearly-MIROC5-rcp85-2015-2100_05000m.nc']);
% lookup table 
lookup = ncload(['../Data/lookup/TaSMB_trans_lookup_b25_MARv3.9-MIROC5-rcp85.nc']);
modscen='M39_MIROC5-rcp85';

% dummy lookup for zero
dummy0 = lookup.aSMB_ltbl(:,1,1);


% Load a geometry for reconstruction
amod = 'OBS';
nc=ncload(['../Models/' amod '/lithk_05000m.nc']);
nc1=ncload(['../Models/' amod '/topg_05000m.nc']);
nc2=ncload(['../Models/' amod '/sftgif_05000m.nc']);
nc.lithk = double(nc.lithk);

% Operate on ice thickness
thi = nc.lithk.*nc2.sftgif;
thi_re = nc.lithk.*nc2.sftgif;
bed = nc1.topg;
sur0 = bed+thi;
sur = bed+thi;
sur_re=bed+thi_re;

mask=nc2.sftgif;

nt=length(lookup.time);

bint_re=zeros(nb,nt);

% output array
tdsmb_re=zeros(size(d0.aSMB));
tsur=zeros(size(d0.aSMB));
tsur_re=zeros(size(d0.aSMB));

%for t=1:5 % year loop
%for t=(nt-5):nt % year loop
for t=1:nt % year loop

    t
    dsd=d0.aSMB(:,:,t) * secpyear / 1000 * rhof/rhoi;
    dsd_re=zeros(size(dsd));

    %% loop through basins
    for b=1:nb
        %% set current basin and lookup
        % Determine surface to calculate remapping on
        if(flg_t==0)
            %% sur0: constant sur
            eval(['sur_b=sur0.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
            eval(['ima_b=mask.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        elseif(flg_t==1)
            %% sur: sur changing with MAR DSMB
            eval(['sur_b=sur.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
            eval(['ima_b=mask.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        elseif(flg_t==2)
            %% sur_re: sur changing with remapped MAR DSMB
            eval(['sur_b=sur_re.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
            eval(['ima_b=mask.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        elseif(flg_t==3)
            %% sur0: constant sur + dSMB/dz (below)
            eval(['sur_b=sur0.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
            eval(['ima_b=mask.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        elseif(flg_t==4)
            %% sur: sur changing with MAR DSMB + dSMB/dz (below)
            eval(['sur_b=sur.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
            eval(['ima_b=mask.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        elseif(flg_t==5)
            %% sur_re: sur changing with remapped MAR DSMB + dSMB/dz (below)
            eval(['sur_b=sur_re.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
            eval(['ima_b=mask.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        else
            disp(['Warning: not an option for flg_t=', num2str(flg_t) ])
            return;
        end
        look0=lookup.aSMB_ltbl(:,wbas.n0(b),t);
         %% set neighbor basin and lookup
        look0 = dummy0;
        if (wbas.n0(b)>0)
            look0=lookup.aSMB_ltbl(:,wbas.n0(b),t);
        end
        look1 = dummy0;
        if (wbas.n1(b)>0)
            look1=lookup.aSMB_ltbl(:,wbas.n1(b),t);
        end
        look2 = dummy0;
        if (wbas.n2(b)>0)
            look2=lookup.aSMB_ltbl(:,wbas.n2(b),t);
        end
        look3 = dummy0;
        if (wbas.n3(b)>0)
            look3=lookup.aSMB_ltbl(:,wbas.n3(b),t);
        end
        look4 = dummy0;
        if (wbas.n4(b)>0)
            look4=lookup.aSMB_ltbl(:,wbas.n4(b),t);
        end
        look5 = dummy0;
        if (wbas.n5(b)>0)
            look5=lookup.aSMB_ltbl(:,wbas.n5(b),t);
        end
        look6 = dummy0;
        if (wbas.n6(b)>0)
            look6=lookup.aSMB_ltbl(:,wbas.n6(b),t);
        end
        
        %% use lookup table to determine DSMB
        dsd_b0 = interp1(lookup.z,look0(:),sur_b);
        dsd_b1 = interp1(lookup.z,look1(:),sur_b);
        dsd_b2 = interp1(lookup.z,look2(:),sur_b);
        dsd_b3 = interp1(lookup.z,look3(:),sur_b);
        dsd_b4 = interp1(lookup.z,look4(:),sur_b);
        dsd_b5 = interp1(lookup.z,look5(:),sur_b);
        dsd_b6 = interp1(lookup.z,look6(:),sur_b);

        if (flg_weigh == 0)
            %% combine according to weights
            dsd_b = dsd_b0.*wbas.wg;
        else
            dsd_b = dsd_b0.*wbas.wgc0 + dsd_b1.*wbas.wgc1 + dsd_b2.*wbas.wgc2 + dsd_b3.*wbas.wgc3 + dsd_b4.*wbas.wgc4 + dsd_b5.*wbas.wgc5 + dsd_b6.*wbas.wgc6;
        end
%    shade(dsd_b)

        %% extended integral dsmb for this basin
        dsd_ex = dsd.*ima_b;
        bint_ex(b)=nansum(nansum(dsd_ex.*af2))*dx*dy;
        
        %% mask
        dsd_b = dsd_b.*mask;

        %% integral dsmb for this basin
        bint_re(b)=nansum(nansum(dsd_b.*af2))*dx*dy;

        if (flg_match == 1)        
            %% adjust to match total
            dsd_b = dsd_b * (lookup.bint(b)/bint_re(b));
        end
        
        %% check integral again 
        bint_out(b)=nansum(nansum(dsd_b.*af2))*dx*dy;

        %% replace nan by zeros to add all basins together
        dsd_b(isnan(dsd_b))=0;
        dsd_re = dsd_re+dsd_b;

    end
    %% end basin loop

    if(flg_t==3 | flg_t==4 | flg_t==5 )
        %% dSMB/dz
        dsmb = d1.dSMBdz(:,:,t) .* (sur-sur0) * rhof/rhoi;
        dsmb_re = d1.dSMBdz(:,:,t) .* (sur_re-sur0) * rhof/rhoi;
        dsd = dsd + dsmb;
        dsd_re = dsd_re + dsmb_re;
    end

    dsd_re = dsd_re * secpyear / 1000 * rhof/rhoi;
    
    %% update surface elevation
    thi_re = max(thi_re + dsd_re, 0);
    sur_re = bed + thi_re;
    thi = max(thi + dsd, 0);
    sur = bed + thi;
    tsur_re(:,:,t) = sur_re;
    tsur(:,:,t) = sur;

    %% collect results
    tdsmb_re(:,:,t) = dsd_re;    

    if (flg_plot) 
        shade_bg(dsd_re)
        colormap(cmap)
        caxis([-4,1])
        print('-dpng', '-r300', ['dsmb_' modscen '_re' sprintf('%02d',t)]) 
        close
        shade_bg(dsd.*mask)
        colormap(cmap)
        caxis([-4,1])
        print('-dpng', '-r300', ['dsmb_' modscen '_or' sprintf('%02d',t)]) 
        close
    end


end
%% end time loop

if(flg_t==0)
    save(['../Data/output/transient0_' modscen, '_' amod],'tsur_re','tsur')
elseif(flg_t==1)
    save(['../Data/output/transient1_' modscen, '_' amod],'tsur_re','tsur')
elseif(flg_t==2)
    save(['../Data/output/transient2_' modscen, '_' amod],'tsur_re','tsur')
elseif(flg_t==3)
    save(['../Data/output/transient3_' modscen, '_' amod],'tsur_re','tsur')
elseif(flg_t==4)
    save(['../Data/output/transient4_' modscen, '_' amod],'tsur_re','tsur')
elseif(flg_t==5)
    save(['../Data/output/transient5_' modscen, '_' amod],'tsur_re','tsur')
end
