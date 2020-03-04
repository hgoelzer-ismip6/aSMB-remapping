# Make SMB anomaly forcing from original MAR files 

set -x

mkdir -p proc
cd proc

amodscen=MIROC5-rcp85
#amodscen=NorESM1-rcp85
#amodscen=MIROC5-rcp26
#amodscen=CSIRO-Mk3.6-rcp85
#amodscen=HadGEM2-ES-rcp85
#amodscen=IPSL-CM5-MR-rcp85
#amodscen=ACCESS1.3-rcp85
#
#amodscen=CESM2-ssp585
#amodscen=CNRM-CM6-ssp126
#amodscen=CNRM-CM6-ssp585
#amodscen=CNRM-ESM2-ssp585
#amodscen=UKESM1-CM6-ssp585

apath=../../Data/dSMB/${amodscen}/aSMB


# Collect files from hist and rcp 
for i in `seq 2091 2100`; do
	cp ${apath}/aSMB_MARv3.9-yearly-${amodscen}-${i}.nc ./MAR_${i}.nc
done

## add time information
for i in `seq 2091 2100`; do
    filename=MAR_${i}.nc
    # yearly sum
    ncks -O --mk_rec_dmn time $filename $filename
    ncap2 -O -s "time=time*0+$i" $filename $filename
done

# concat them to one time series
ncrcat -O  MAR_2091.nc MAR_2092.nc MAR_2093.nc MAR_2094.nc MAR_2095.nc MAR_2096.nc MAR_2097.nc MAR_2098.nc MAR_2099.nc MAR_2100.nc aSMB_MARv3.9-yearly-${amodscen}-2091-2100.nc

# Long-term average
ncra -O  aSMB_MARv3.9-yearly-${amodscen}-2091-2100.nc aSMB_MARv3.9-yearly-${amodscen}_ltm2091-2100.nc
# remove MAPPING artefact
ncks -C -O -x -v MAPPING aSMB_MARv3.9-yearly-${amodscen}_ltm2091-2100.nc aSMB_MARv3.9-yearly-${amodscen}_ltm2091-2100.nc

# regrid
# Path to GDFs
gdfs=../../Data/Grid
ingdf=grid_ISMIP6_GrIS_01000m.nc
outgdf=grid_ISMIP6_GrIS_05000m.nc
cdo remapycon,${gdfs}/${outgdf} -setgrid,${gdfs}/${ingdf} aSMB_MARv3.9-yearly-${amodscen}_ltm2091-2100.nc aSMB_MARv3.9-yearly-${amodscen}_ltm2091-2100_e05000m.nc

# Copy to destination
/bin/cp aSMB_MARv3.9-yearly-${amodscen}_ltm2091-2100_e05000m.nc ../../Data/RCM/

# Clean up
/bin/rm MAR_*.nc
