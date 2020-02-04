# aSMB-remapping
Matlab/Shell workflow for remapping Greenland SMB anomalies

# Prepare input data 

### Copy input files from external Data archive 
```cp <ExtArchive>/ISMIP6_Extensions_05000m.nc ../Data/Basins/```

```cp <ExtArchive>/af2_ISMIP6_GrIS_05000m.nc ../Data/Grid/```

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


### Prepare MAR data (done once per scenario) 

`./process_MAR_reference.sh`

`matlab`

% Build forcing time slice

`save_DSMB.m`

% Build a forcing time series 

`save_trans_DSMB.m`

