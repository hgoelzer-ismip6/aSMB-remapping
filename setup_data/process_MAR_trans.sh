# Make transient SMB anomaly forcing 

set -x

mkdir -p proc
cd proc

avar=aSMB 
#avar=dSMBdz 

amodscen=MIROC5-rcp85

apath=../../Data/dSMB/${amodscen}/${avar}


# Collect files 
for i in `seq 2015 2100`; do
	cp ${apath}/${avar}_MARv3.9-yearly-${amodscen}-${i}.nc ./MAR_${i}.nc
done

## add time information
for i in `seq 2015 2100`; do
    filename=MAR_${i}.nc
    # yearly sum
    ncks -O --mk_rec_dmn time $filename $filename
    ncap2 -O -s "time=time*0+$i" $filename $filename
done

# concat them to one time series
ncrcat -O MAR_2015.nc MAR_2016.nc MAR_2017.nc MAR_2018.nc MAR_2019.nc MAR_2020.nc MAR_2021.nc MAR_2022.nc MAR_2023.nc MAR_2024.nc MAR_2025.nc MAR_2026.nc MAR_2027.nc MAR_2028.nc MAR_2029.nc MAR_2030.nc MAR_2031.nc MAR_2032.nc MAR_2033.nc MAR_2034.nc MAR_2035.nc MAR_2036.nc MAR_2037.nc MAR_2038.nc MAR_2039.nc MAR_2040.nc MAR_2041.nc MAR_2042.nc MAR_2043.nc MAR_2044.nc MAR_2045.nc MAR_2046.nc MAR_2047.nc MAR_2048.nc MAR_2049.nc MAR_2050.nc MAR_2051.nc MAR_2052.nc MAR_2053.nc MAR_2054.nc MAR_2055.nc MAR_2056.nc MAR_2057.nc MAR_2058.nc MAR_2059.nc MAR_2060.nc MAR_2061.nc MAR_2062.nc MAR_2063.nc MAR_2064.nc MAR_2065.nc MAR_2066.nc MAR_2067.nc MAR_2068.nc MAR_2069.nc MAR_2070.nc MAR_2071.nc MAR_2072.nc MAR_2073.nc MAR_2074.nc MAR_2075.nc MAR_2076.nc MAR_2077.nc MAR_2078.nc MAR_2079.nc MAR_2080.nc MAR_2081.nc MAR_2082.nc MAR_2083.nc MAR_2084.nc MAR_2085.nc MAR_2086.nc MAR_2087.nc MAR_2088.nc MAR_2089.nc MAR_2090.nc MAR_2091.nc MAR_2092.nc MAR_2093.nc MAR_2094.nc MAR_2095.nc MAR_2096.nc MAR_2097.nc MAR_2098.nc MAR_2099.nc MAR_2100.nc ${avar}_MARv3.9-yearly-${amodscen}-2015-2100.nc

# Regrid using GDFs
gdfs=../../Data/Grid
ingdf=grid_ISMIP6_GrIS_01000m.nc
outgdf=grid_ISMIP6_GrIS_05000m.nc
cdo remapycon,${gdfs}/${outgdf} -setgrid,${gdfs}/${ingdf} ${avar}_MARv3.9-yearly-${amodscen}-2015-2100.nc  ${avar}_MARv3.9-yearly-${amodscen}-2015-2100_05000m.nc

# Copy to destination
/bin/cp ${avar}_MARv3.9-yearly-${amodscen}-2015-2100_05000m.nc ../../Data/RCM/

# Clean up
/bin/rm MAR_*.nc

