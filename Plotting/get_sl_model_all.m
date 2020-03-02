% get SL contribution from all initMIP models

clear

%modscen='MIROC5_rcp85';
%modscen='initMIP';
modscen='M39_MIROC5_rcp85';

% area factors
da = ncload('../Data/Grid/af2_ISMIP6_GrIS_05000m.nc');
af2 = double(da.af2(:,:));

%% dim
dx=5000;dy=5000;

%% model individual parameters
load ../Data/initMIP/ch_A5.mat

%% output
%iadsmb_re=zeros([100,ch.n]);
iadsmb_re=zeros([86,ch.n]);

for m=1:ch.n
    %for m=17:19

    amod = ch.ids{ch.order(m)};

    rhoi=ch.rhoi{1};
    v2mmsl=-1e-12/361.8*rhoi;
    scl=v2mmsl;

    %% models
    t0=load(['../Data/initMIP/transient0_' modscen '_' amod '.mat']);
    % remapped
    tsur_ref = t0.tsur_re(:,:,1);
    tsur_ref3= repmat(tsur_ref,[1,1,size(t0.tsur_re,3)]);
    tdsur = t0.tsur_re-tsur_ref3;
    iadsmb_re(:,m)=squeeze(nansum(nansum(tdsur.*af2))*dx*dy*scl);
    % original
    tsur_ref = t0.tsur(:,:,1);
    tsur_ref3= repmat(tsur_ref,[1,1,size(t0.tsur,3)]);
    tdsur = t0.tsur-tsur_ref3;
    iadsmb(:,m)=squeeze(nansum(nansum(tdsur.*af2))*dx*dy*scl);

end

% update
prog.n = ch.n;
prog.igrp = ch.igrp;
prog.imod = ch.imod;
prog.igrpmod = ch.igrpmod;
prog.iadsmb_re = iadsmb_re;
prog.iadsmb = iadsmb;
prog.ids = ch.ids;

% save
save(['../Data/initMIP/A5_iadsmb_trans_' modscen], 'prog');
