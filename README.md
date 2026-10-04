# Black Hole Imaging and Ray Tracing

Numerical astrophysics project developed as part of the MSc in Astrophysics,
Space Science and Planetology (ASEP) at the University of Toulouse.

**Authors:** 
[@Willydrogene](https://github.com/Willydrogene) and
[@Zoltrak-Kiruwa](https://github.com/Zoltrak-Kiruwa)

## Project Overview

This project explores the numerical simulation of photon trajectories around
a Schwarzschild black hole using relativistic geodesics and ray tracing.

The numerical integration was implemented in C using a fourth-order
Runge-Kutta (RK4) method, while Python was used for data processing and
visualisation.

The project includes:

- numerical integration of photon trajectories around a Schwarzschild black hole;
- ray tracing and reconstruction of the observed black hole image;
- modelling of an accretion disk;
- gravitational and Doppler redshift effects;
- higher-order images produced by strongly deflected photon trajectories;
- reconstruction of spectral line profiles for different disk inclinations;
- CPU parallelisation using OpenMP;
- an additional gravitational lensing simulation.

## Report

The complete report for the main black hole imaging project is available here:

**[Black_Hole_Report.pdf](./Black_Hole_Report.pdf)**

The gravitational lensing extension was developed after the main project and
is therefore **not included in the report**. However, the corresponding C and
Python source files are available in this repository.

## Gravitational Lensing Extension

In addition to the original project, we implemented a gravitational lensing
simulation allowing an arbitrary background image to be distorted by the
black hole.

The source files for this extension can be found in this repository
(`lentille.c`, `lentille_xy.c`, and the associated Python visualisation scripts).

Example images and simulation outputs are available here:

**[Gravitational lensing images and outputs](https://drive.google.com/drive/folders/1xOz_J9GVHaxCri4efSZySkAkuWe9v_0y?usp=sharing)**

## Technologies

- **C** — numerical integration and ray tracing
- **Python** — data processing and visualisation
- **OpenMP** — CPU parallelisation and performance optimisation
- **RK4** — numerical integration of relativistic geodesics

## Supervisors

Laurène Jouve and Douglas Marshall  
Université de Toulouse
