# TractionModel2D

Model described in: Oakes, P. W., Banerjee, S., Marchetti, M. C., & Gardel, M. L. (2014). Geometry regulates traction stresses in adherent cells. Biophysical journal, 107(4), 825-833. Developed by Shiladitya Banerjee.

This repository contains finite element simulation scripts, written in Matlab, for the cell tracton model described in Oakes et al Biophys J 2014. Simulations require installation of the Partial Differential Equation Toolbox in Matlab.

The repository contains four directories that provide simulation scripts for traction stress modeling for (1) cells in circular micropatterns, (2) stadium-shaped micropatterns, (3) Fibroblast cell shapes and (4) Any input cell shapes.

Model parameters are specific to each cell type and substrate properties, and must be calibrated with experimental data to yield quantitatively accurate data.

To simulate the traction stress map, simply run the script 'pde*.m' in each directory. The script is commented in 'any_shape/pdesol_uc.m'.

Model Parameters:

E - Young's modulus of the cell

```math \nu ```- Poisson' ratio of the cell
