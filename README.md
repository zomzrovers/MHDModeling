# MHD Modeling of Phase Change in PCMs

**Mentor:** Prof. Virkeshwar Kumar, Department of Mechanical Engineering, IIT Kanpur  
**Duration:** May 2025 – Present

## Overview

This project investigates the magnetohydrodynamic (MHD) behavior of phase change materials (PCMs) under the influence of magnetic fields during melting and solidification. By modeling the coupled phenomena of fluid flow, heat transfer, and phase change, the project aims to understand how magnetic fields can be used to optimize thermal control in advanced manufacturing and thermal energy storage systems.

## Objective

- Quantify the impact of magnetic field strength and orientation on melting and solidification dynamics in PCMs.
- Generate insights to support thermal management strategies in real-world applications using materials like gallium and paraffin wax.

## Methodology

- Developed a custom MHD solver in OpenFOAM by extending phase change and buoyancy-driven flow solvers to include electromagnetic body forces.
- Conducted over 50 simulations varying magnetic field strength and orientation.
- Studied the effects on flow structure, solid fraction evolution, and dominant heat transfer mode (conduction vs. convection).
- Post-processed data to visualize temperature fields, velocity magnitudes, and melting front progression.

## Key Findings

- Optimal magnetic field orientation led to a 31% and 17% reduction in solid fraction during melting and solidification, respectively.
- Magnetic fields suppressed natural convection, reducing fluid velocity and promoting conduction-dominated heat transfer.

## Simulation Setup

- **Software**: OpenFOAM (custom solvers and boundary conditions)
- **Geometry**: 2D rectangular and axisymmetric domains
- **Materials**: Gallium and Paraffin Wax
- **Boundary Conditions**: Prescribed heat flux / temperature gradient; applied magnetic fields with different orientations
- **Mesh Resolution**: High-resolution grids (up to 100k+ cells) for capturing phase interfaces accurately


