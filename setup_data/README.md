# aSMB-remapping
Matlab/Shell workflow for remapping Greenland SMB anomalies

# Prepare input data 

### Copy input files from external Data archive 
```cp <ExtArchive>/ISMIP6_Extensions_05000m.nc ../Data/Basins/```

```cp <ExtArchive>/af2_ISMIP6_GrIS_05000m.nc ../Data/Grid/```

```cp <ExtArchive>/zmask_05000m.nc ../Data/RCM/```

```cp <ExtArchive>/orog_05000m.nc ../Models/OBS/```

```cp <ExtArchive>/sftgif_05000m.nc ../Models/OBS/```

```cp <ExtArchive>/orog_05000m.nc ../Models/VUBGISM/```

```cp <ExtArchive>/sftgif_05000m.nc ../Models/VUBGISM/```

### Prepare Basins (done only once)

`./make_masks.sh`

`matlab`

% Save basins in useful mask format 

`save_extbasins.m`

% Calculate basin weights

`save_extbasin_neighbour.m`

`save_extbasin_scale_div.m`


### Prepare RCM data (done once per scenario) 

% Make SMB reference

`./process_MAR_reference_01.sh`

% Calcualte SMB anomalies and gradients

`matlab`
`save_trans_DSMB_01.m`

% Build a forcing time series 

`./process_MAR_trans.sh`

% Build forcing time slice

`./process_MAR.sh`

