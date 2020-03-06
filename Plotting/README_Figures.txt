README
Figure numbers for revised manuscript version 035

Most scripts and material in
/Volumes/ISMIP6/Projections/Greenland/SMB_ISMIP6_25


Figure 01
ELAscheme.pptx


Figure 02
% Define basins
save_extbasins
% Plot basins 
plot_wbas_div_bw


Figure 03
dsmb_5km_basin_fit
maximise to window, then print
print -dpng -r300 bfit_gradients
-->bfit_gradients_nb25_MARv37_MIROC5_rcp85.png


Figure 04
% Define basins
save_extbasins
% Plot basins 
plot_wbas_div_bw

Figure 05
reconstruct_aSMB
iism = 0;
-->
ddsmb_div_M39_MIROC5_rcp85_map-ext_obs_wgt1.png
dsmb_div_M39_MIROC5_rcp85_mapmod_obs_wgt1.png
dsmb_div_M39_MIROC5_rcp85_orgobs.png
Then Figure5_Set.pptx
save as png
--> Figure5_Set.png

Figure 06
reconstruct_aSMB

iism = 0;
--> dsmb_basinint_M39_MIROC5_rcp85_re.png

Figure 07
reconstruct_aSMB
iism = 1; VUB
-->
ddsmb_div_M39_MIROC5_rcp85_ext-org_vubgism_wgt1.png
ddsmb_div_M39_MIROC5_rcp85_map-ext_vubgism_wgt1.png
dsmb_div_M39_MIROC5_rcp85_mapmod_vubgism_wgt1.png
dsmb_div_M39_MIROC5_rcp85_orgmod_vubgism.png
dsmb_div_M39_MIROC5_rcp85_orgobs.png
dsur_divM39_MIROC5_rcp85_re.png
Then Figure07_Set.pptx
save as png
--> Figure07_Set.png


Figure 08
reconstruct_aSMB
iism = 1; VUB
-->
dsmb_basinint_M39_MIROC5_rcp85mod_map_vubgism_wgt1.png
ddsmb_basinint_M39_MIROC5_rcp85mod_map_vubgism_wgt1.png


Figure 09
plot_trans_lookup
maximise than manually print to png
--> trans_lookup_MAR39_MIROC5_rcp85_b25.png
crop margins in preview


Figure 10
meta_trans_dsmb
calls trans_dsmb_5km_basin_reconstruct_dsdz
plot_trans_results_dsdz
-->
Dh_c_2100_MAR39_obs.png
Dh_0_2100_MAR39_obs.png
Dh_param_MAR39_obs.png
Then Figure10_Set.pptx
Print, select paper screen19x9
save as png
--> Figure10_Set.png



Figure 11
meta_trans_dsmb
calls trans_dsmb_5km_basin_reconstruct_dsdz
plot_trans_results_dsdz
-->
Dh_2_2100_MAR39_obs.png
Dh_lapse_2100_MAR39_obs.png
Dh_lapseD_2100_MAR39_obs.png
Then Figure11_Sets.pptx
Print, select paper screen19x9
save as png
--> Figure11_Sets copy.png


Figure 12
trans_dsmb_5km_basin_reconstruct_all
get_sl_model_all
bar_sl_model_all
-->
A5_bar_trans_M39_MIROC5_rcp85.png
A5_bar_trans_remap_M39_MIROC5_rcp85.png


=========== Supplement ===============

Figure S3
reconstruct_aSMB
iism = 2; iism = 3; iism = 4;
-->
dsmb_div_M39_MIROC5-rcp85_mapmod_bgcbi_wgt1.png
dsmb_div_M39_MIROC5-rcp85_mapmod_jplissm_wgt1.png
dsmb_div_M39_MIROC5-rcp85_mapmod_mpimpism_wgt1.png

# Figure S4
reconstruct_aSMB
iism = 2; iism = 3; iism = 4;
-->
dsmb_basinint_M39_MIROC5-rcp85mod_map_bgcbi_wgt1.png
dsmb_basinint_M39_MIROC5-rcp85mod_map_jplissm_wgt1.png
dsmb_basinint_M39_MIROC5-rcp85mod_map_mpimpism_wgt1.png
