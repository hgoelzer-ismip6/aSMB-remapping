# aSMB-remapping
Matlab/Shell workflow for remapping Greenland SMB anomalies

# Reconstruct DSMB using lookup table

# Workflow
### Use time dependent lookup table to reconstruct aSMB and dSMBdz at model surface

`matlab`

`reconstruct_aSMB.m`

`trans_reconstruct_aSMB.m`

`trans_reconstruct_dSMBdz.m`

--> ../Models/<MODEL>/<scenario>/aSMB/aSMB_MARv3.9-yearly-<scenario>-<MODEL>-nnnn.nc

--> ../Models/<MODEL>/<scenario>/dSMBdz/dSMBdz_MARv3.9-yearly-<scenario>-<MODEL>-nnnn.nc


### Run schematic scenario
`meta_trans_dsmb.m`

`trans_dsmb_5km_basin_reconstruct_dsdz.m`
