# Seis3D

Compute one way coupled radial pressure diffusion in poroelastic medium, 3D aseismic fault slip, and off-fault seismicity

Taiyi A. Wang, Eric M. Dunham, 2022.

## Usage:
 
projects/
 - cooperBasin.m: computes all required information in the paper: "Hindcasting injection-induced aseismic slip and microseismicity at the Cooper Basin Enhanced Geothermal Systems Project"
 - mk_plots: make all plots in the paper

source_code/
 - functions that scripts in /projects will call
 
benchmark/
 - benchmark static elasticity solver and pressure diffusion solver

Note: 
1. All input parameter values are set up in setup_model.m under /source_code
2. All input data (injection volume rate etc.) are setup in setup_data.m under /source_code

Warning:
cooperBasin.m may take a few days to run to complete. The most time consuming step is computing aseismic slip on main fault.
