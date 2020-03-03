# aSMB-remapping
Matlab/Shell workflow for remapping Greenland SMB anomalies

# Make SMB lookup table

### Create time slice lookup table 

`basin_fit_aSMB.m`
--> Data/lookup/aSMB_lookup_b25_MARv3.9-<scenario>.nc

### Create time dependent lookup table for a given RCM simulation

`trans_basin_fit_aSMB.m`
--> Data/lookup/TaSMB_trans_lookup_b25_MARv3.9-<scenario>.nc

`trans_basin_fit_dSMBdz.m`
--> Data/lookup/TdSMBdz_trans_lookup_b25_MARv3.9-<scenario>.nc

