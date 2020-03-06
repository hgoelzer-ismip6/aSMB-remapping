% Fit for a number of products and basins

addpath('../toolbox')

% Plot settings
colororder = distinguishable_colors(100);
set(0,'DefaultAxesColorOrder', colororder)
set(0,'DefaultTextInterpreter','none');
set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFIgurePaperPositionMode','auto')
set(0, 'DefaultAxesFontSize', 16)
set(0, 'DefaultLineLineWidth', 2)

% comment clear, mul and iscen in script dsmb_5km_basin_fit20xN.m
for iscen = 1:6
    for mul =1:5
        clear lookup
        aSMB_5km_basin_fit20xN
        close 
    end
end
