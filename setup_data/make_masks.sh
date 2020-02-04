#!/bin/bash
# Save a series of basin masks

for id in `seq 1 25`; do
    echo ${id}
    ncap2 -A -s "basin${id}=IDs*0.; where(IDs==${id}) basin${id}=1" -v ../Data/Basins/ISMIP6_Extensions_05000m.nc ../Data/Basins/ISMIP6Masks25_05000m.nc
done
