# aSMB-remapping
Matlab/Shell workflow for remapping Greenland SMB anomalies

The input dataset needed to reproduce the results (<ExtArchive>) can be found on Zenodo: doi:10.5281/zenodo.3760526

# Prepare input data 

### Copy input files from external Data archive 

`cp <ExtArchive>/ISMIP6_Extensions_05000m.nc ../Data/Basins/`

`cp <ExtArchive>/af2_ISMIP6_GrIS_05000m.nc ../Data/Grid/`
`cp <ExtArchive>/days_1900-2300.txt ../Data/Grid/`
`cp <ExtArchive>/grid_ISMIP6_GrIS_01000m.nc ../Data/Grid/`
`cp <ExtArchive>/grid_ISMIP6_GrIS_05000m.nc ../Data/Grid/`
`cp <ExtArchive>/zmask_05000m.nc ../Data/Grid/`

`cp <ExtArchive>/grid_MARv3.9_05000m.nc ../Data/RCM/`
`cp <ExtArchive>/orog_MARv3.9_05000m.nc ../Data/RCM/`

`cp <ExtArchive>/OBS/orog_05000m.nc ../Models/OBS/`
`cp <ExtArchive>/OBS/sftgif_05000m.nc ../Models/OBS/`
`cp <ExtArchive>/OBS/lithk_05000m.nc ../Models/OBS/`
`cp <ExtArchive>/OBS/topg_05000m.nc ../Models/OBS/`

`cp <ExtArchive>/initMIP/VUBGISM/orog_05000m.nc ../Models/VUBGISM/`
`cp <ExtArchive>/initMIP/VUBGISM/sftgif_05000m.nc ../Models/VUBGISM/`

`cp <ExtArchive>/initMIP/fi_A5.mat ../Data/initMIP`
`cp <ExtArchive>/initMIP/ch_A5.mat ../Data/initMIP`


### Prepare Basins (done only once)

`./make_masks.sh`

`matlab`

% Save basins in mask format 

`save_extbasins.m`

% Calculate basin weights

`save_extbasin_neighbour.m`

`save_extbasin_scale_div.m`

--> Results in ../Data/Basins

### Prepare RCM data (done once per scenario) 

The MAR data originally comes as yearly SMB and dRUNdz at 1 km resolution. E.g.
ftp://ftp.climato.be/fettweis/MARv3.9/ISMIP6/GrIS/MIROC5-histo_1950_2005
We calculate anomalies against the 1960-1989 reference period and interpolate to 
5 km resolution. The result of this operation (the 5 km files) are provided for 
convenience. In that case the scripts below are not needed.

% Make SMB reference (~10 min)

`./process_MAR_reference_01.sh`

% Calcualte SMB anomalies and gradients (~60 min)

`matlab`
`save_trans_DSMB_01.m`

% Build a forcing time series (2x ~10 min)
% Run twice: for avar=aSMB and avar=dSMBdz

`./process_MAR_trans.sh`

% Build forcing time slice (~1 min)

`./process_MAR.sh`

--> Results in ../Data/RCM
