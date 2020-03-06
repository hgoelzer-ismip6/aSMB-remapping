% reconstruct for a number of products and basins

% comment clear,  mul and iscen in script dsmb_5km_basin_reconstruct_scale_div20xN
rmss = zeros(6,5);
for iscen = 1:6
    for mul =1:5
        aSMB_5km_basin_reconstruct_scale_div20xN
        %% remember error output
        rmss(iscen,mul) = rmsn;
    end
end

save meta_rms_M39 rmss
