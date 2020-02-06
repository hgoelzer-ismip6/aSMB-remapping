% Extract aSMB data

clear

%%% MAR 3.9
%pdata='/Volumes/Storage/ISMIP6_Disk/Data/GrIS/MAR/MAR3.9/MAR_05000m/';
pdata='/Volumes/ISMIP6/Data/GrIS/MAR/MAR3.9/MAR_05000m/';

d0=ncload([pdata 'MARv3.9-yearly-MIROC5-ref_05000m.nc']);  
d1=ncload([pdata 'MARv3.9-yearly-MIROC5-2015-2100_05000m.nc']);  
DSMB=(mean(d1.SMB(:,:,77:86),3)-d0.SMB(:,:))/1000;
% other info
LAT=d0.LAT;
MSK=d0.MSK;
SH=d0.SRF;
% write out
save ../Data/RCM/DSMB_MARv3.9_MIROC5_rcp85.mat DSMB LAT MSK SH
