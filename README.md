# SiphonInvagination_VertexModel
MATLAB code for vertex model simulations of Ciona atrial siphon invagination. Paper: [Bidirectional redistribution of actomyosin drives epithelial invagination in ascidian siphon tube morphogenesis](https://doi.org/10.7554/eLife.108588.2), J. Qiao, P. Yu, H. Peng, W. Shi, B. Li, and B. Dong, eLife (2026).

## Overview

This repository contains the MATLAB implementation of the two-dimensional vertex model used to investigate the mechanical mechanisms underlying Ciona atrial siphon invagination. The model represents the epithelial tissue as polygonal cells arranged in a circular geometry and incorporates cell area elasticity, passive cortical contractility, active apical/basal/lateral tensions, and tissue bending elasticity.

The default parameter set in `Main.m` corresponds to the representative simulation reported in the manuscript. Other simulation conditions were generated using the same model framework by modifying the corresponding parameters according to the values described in the manuscript.

## Code organization

- `Main.m` — Main simulation script. It defines model parameters, prescribed myosin dynamics, active tensions, time integration, and simulation output.

- `Mesh.m` — Initializes the epithelial geometry and vertex network, including cell connectivity and identification of apical, basal, and lateral edges.

- `ForceGroup.m` — Calculates the mechanical forces acting on each vertex, including contributions from cell area elasticity, passive cortical contraction, active tensions, lumen area constraint, and tissue bending elasticity.

- `CellSize.m` — Calculates cell geometrical properties, including cell area, perimeter, lumen area, and cell-center positions.

- `SaveData.m` — Stores simulation configurations and model parameters during the simulation.

## Requirements

The code was developed in MATLAB R2021a and requires the Statistics and Machine Learning Toolbox.

Before running the simulation, place all `.m` files in the same directory.

## Copyright and Contact

This code was developed by Dr. Pengyu Yu, from Bo Li Lab, Department of Engineering Mechanics, Tsinghua University.

Copyright © 2026 Bo Li Lab, Tsinghua University.

For questions regarding the associated publication, please contact the corresponding authors listed in the paper.
