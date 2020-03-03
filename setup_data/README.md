# aSMB-remapping
Matlab/Shell workflow for remapping Greenland SMB anomalies

# Prepare input data 

### Copy input files from external Data archive 
```cp <ExtArchive>/ISMIP6_Extensions_05000m.nc ../Data/Basins/```

```cp <ExtArchive>/af2_ISMIP6_GrIS_05000m.nc ../Data/Grid/```

```cp <ExtArchive>/zmask_05000m.nc ../Data/RCM/```

```cp <ExtArchive>/OBS/orog_05000m.nc ../Models/OBS/```

```cp <ExtArchive>/OBS/sftgif_05000m.nc ../Models/OBS/```

```cp <ExtArchive>/OBS/lithk_05000m.nc ../Models/OBS/```

```cp <ExtArchive>/OBS/topg_05000m.nc ../Models/OBS/```

```cp <ExtArchive>/VUBGISM/orog_05000m.nc ../Models/VUBGISM/```

```cp <ExtArchive>/VUBGISM/sftgif_05000m.nc ../Models/VUBGISM/```


### Prepare Basins (done only once)

`./make_masks.sh`

`matlab`

% Save basins in mask format 

`save_extbasins.m`

% Calculate basin weights

`save_extbasin_neighbour.m`

`save_extbasin_scale_div.m`


### Prepare RCM data (done once per scenario) 

The MAR data originally comes as yearly SMB and dRUNdz at 1 km resolution. 
We calculate anomalies against the 1960-1989 reference period and interpolate to 
5 km resolution. The result of this operation (the 5 km files) are provided for 
convenience. In that case the scripts below are not needed.

% Make SMB reference

`./process_MAR_reference_01.sh`

% Calcualte SMB anomalies and gradients

`matlab`
`save_trans_DSMB_01.m`

% Build a forcing time series 

`./process_MAR_trans.sh`

% Build forcing time slice

`./process_MAR.sh`

