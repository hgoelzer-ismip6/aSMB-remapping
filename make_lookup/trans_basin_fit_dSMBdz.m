% Create time dependent lookup table for dSMBdz
% For remapping this is based on runoff gradients dRUNdz!

clear

addpath('../toolbox')

%% Settings
rcm = 'MARv3.9';

gcm = 'MIROC5';
scen = 'rcp85';

%%%%%%%

% replace nan for high elevation with last finite value?
flg_nanfill = 1;

colors=get(0,'DefaultAxesColorOrder');

% basin definition
load ../Data/Basins/ExtBasinMasks25_05000m.mat
nb=length(bas.ids);

% area factors
da = ncload('../Data/Grid/af2_ISMIP6_GrIS_05000m.nc');
af = double(da.af2(:,:));

% res
dx=5000;dy=5000;

inpath = ['../Data/dSMB/' gcm '-' scen ];

% MAR surface 
d0 = ncload(['../Data/RCM/MARv3.9_topg_05000m.nc']);
sur = d0.topg(:,:);

% scenario specific 
lookup_file = ['trans_lookup_b25_MARv3.9-' gcm '-' scen ];

d1 = ncload(['../Data/RCM/dSMBdz_MARv3.9-yearly-' gcm '-' scen '-2015-2100_05000m.nc']);

% timer
time = 2015:2100;
nt = length(time);

ds = ncload('../Models/OBS/sftgif_05000m.nc');
mask = double(ds.sftgif);

%% fit a lookup table to the data
%% centers at sstep intervals with ds range
sstep = 100;
ds = 100;
ss=0:sstep:3500;
ns=length(ss);

table=zeros([nb,ns,nt]);
bint=zeros(nb,nt);

oldlook = zeros(2,length(ss)+1); 

msg = (['running year, basin: 00,00']);
fprintf(msg);
%for t = 1 % year loop
%for t = 1:5 % year loop
%for t = nt % year loop
for t = 1:nt % year loop

%    t
    fprintf(['\b\b\b\b\b']);
    fprintf([sprintf('%02d',t), ',00']);
    dSMBdz = d1.dSMBdz(:,:,t);

    
%    figure
    for b = 1:nb
        
        fprintf(['\b\b\b']);
        fprintf([',' sprintf('%02d',b)]);
        eval(['dSMBdz_b=dSMBdz.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        eval(['sur_b=sur.*(bas.basin' num2str(b) './bas.basin' num2str(b) ');']);
        
%        subplot(5,4,b)
%        hold on; box on;
%        plot(sur_b(:),dSMBdz_b(:),'.');

        %% integral dsmb for this basin
        bint(b,t)=nansum(nansum(dSMBdz_b.*af.*mask))*dx*dy;
        
        %% fit a lookup table to the data
        %% centers at sstep intervals with ds range
        look = zeros(2,ns); 
        n=1;
        yold=0.;
        for s0=ss(2:end)
            n=n+1;
            %% find local average DSMB for given elevation range
%            [y1,ysel,s1,ssel]= find_local_average(sur_b,dSMBdz_b,s0,ds);
            [y1,ysel,s1,ssel]= find_local_median(sur_b,dSMBdz_b,s0,ds);
            %% fill nans with last value
            if(isnan(y1) && flg_nanfill)
                look(:,n)=[s0,yold];
%                plot(s0,yold,'o','Color',colors(n,:))
            else
                look(:,n)=[s0,y1];
                yold=y1;
%                plot(s0,y1,'o','Color',colors(n,:))
            end
%            plot(ssel,ysel,'.','Color',colors(n,:))
        end
        %% fill first value for h=0 from second
        look(2,1) = look(2,2);

        %plot(look(1,:),look(2,:),'-k')
        %axis([0 3300 -5 2])
        %axis([0 3300 -2 1])
        %title(['B' num2str(bas.ids{b}) ' t:' num2str(t)])

        if(flg_nanfill)
            % if first values are zero, set neighbor basin
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
        end
    
        table(b,:,t)=look(2,:);

    end % end basin loop

    %% manual correction for some basins needed?

end % end year loop
fprintf('\n');

%% Write netcdf
nz = length(ss);
td = time;
nt = length(td);

% permute to get dsmb_table(h,b,t)
table_out = permute(table,[2,1,3]);

% Write netcdf
ancfile = ['../Data/lookup/TdSMBdz_' lookup_file '.nc'];
ncwrite_TDSMBTable(ancfile,ss,td,table_out,'dSMBdz_ltbl',{'z','b','time'});
nccreate(ancfile,'bint','Dimensions',{'b',nb,'time',nt}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'bint',bint);
