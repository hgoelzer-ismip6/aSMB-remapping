# aSMB-remapping
Matlab/Shell workflow for remapping Greenland SMB anomalies

# Sensitivito to number of basins

# make lookup table
`meta_dsmb_5km_basin_fit.m`
calls	`dsmb_5km_basin_fit20xN.m`

# reconstruct SMB for a number of basins
`meta_dsmb_5km_basin_reconstruct.m`
calls	`dsmb_5km_basin_reconstruct_scale_div20xN.m`

# Plot error = f(nbasin)
`plot_metaxN.m`
