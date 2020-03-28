function ncwrite_GrIS(ancfile,avar,varname,dimnames,res)
% write variables to a netcdf file and add ISMIP6 specific information
% function ncwrite_GrIS(ancfile,avar,varname, dimnames, res)

% TODO: check input arguments and sizes

% Always overwrite, remove destination if exists
if exist(ancfile, 'file') ~= 0;
    delete(ancfile)
end

% Overwrite coordinate dim names assuming standard ordering (x, y, ...)
dimnames{1}='x1';
dimnames{2}='y1';

% Open netCDF file.
ncid = netcdf.create(ancfile,'CLOBBER');

% Define the dimensions of the variable.
dimids=[];
for i=1:ndims(avar)
  dimids = [dimids,netcdf.defDim(ncid,dimnames{i},size(avar,i))];
end
% Define a new variable in the file.
my_varID = netcdf.defVar(ncid,varname,'double',dimids);

% Leave define mode and enter data mode to write data.
netcdf.endDef(ncid);

% Write data to variable.
netcdf.putVar(ncid,my_varID,avar);

netcdf.close(ncid);

% Standard ISMIP6 coordinates
dx = -713000;
dy = -3396000;
nx1=1664;
ny1=2816;
nx=(nx1)/res;
ny=(ny1)/res;
xd=single(zeros(nx,1));
yd=single(zeros(ny,1));
for ip=1:nx
    xd(ip) = (dx + (ip-1) * res*1000);
end
for jp=1:ny
    yd(jp) = (dy + (jp-1) * res*1000);
end

%% Coordinates
nccreate(ancfile,'x1','Dimensions',{'x1',nx}, 'Datatype','single', 'Format','classic');
nccreate(ancfile,'y1','Dimensions',{'y1',ny}, 'Datatype','single', 'Format','classic');
ncwrite(ancfile,'x1',xd);
ncwrite(ancfile,'y1',yd);
ncwriteatt(ancfile,'x1', 'units', 'm') ;
ncwriteatt(ancfile,'y1', 'units', 'm') ;
ncwriteatt(ancfile,'x1', 'standard_name', 'projection_x_coordinate') ;
ncwriteatt(ancfile,'y1', 'standard_name', 'projection_y_coordinate') ;
ncwriteatt(ancfile,'x1', 'axis', 'x') ;
ncwriteatt(ancfile,'y1', 'axis', 'y') ;

%% Attributes
ncwriteatt(ancfile,'/','proj4','+init=epsg:3413')
ncwriteatt(ancfile,'/','Description',['Created by Heiko Goelzer (h.goelzer@uu.nl), IMAU, ', date ])

% Mapping
ncwrite_mapping_GrIS(ancfile);
