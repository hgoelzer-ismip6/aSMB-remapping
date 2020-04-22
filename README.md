# aSMB-remapping
Matlab/Shell workflow for remapping Greenland SMB anomalies

The input dataset needed to reproduce the results can be found on Zenodo:
doi:10.5281/zenodo.3760526

The remapping approach is documented in the following publication:
Goelzer, H., Noel, B. P. Y., Edwards, T. L., Fettweis, X., Gregory, J. M., Lipscomb, W. H., van de Wal, R. S. W., and van den Broeke, M. R.: Remapping of Greenland ice sheet surface mass balance anomalies for large ensemble sea-level change projections, The Cryosphere Discuss., https://doi.org/10.5194/tc-2019-188, in review, 2019.


# Workflow
### create input data (once and once per scenario) 
`setup_data/`

### Get modelled ice geometry and mask (once per initial state)
`Models/`

### Build lookup table (once per scenario)
`make_lookup/`

### Apply SMB scenario using lookup table 
`aSMB_remap/`



### Plotting scripts
`Plotting/`


# Input data
`Data/`

# matlab toolbox
`toolbox/`

