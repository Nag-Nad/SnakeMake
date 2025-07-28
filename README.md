# MD Simulation Workflow with Snakemake

This repository contains a Snakemake workflow for running molecular dynamics (MD) simulations using GROMACS and post-processing for ensemble docking.

## 🔧 File Description

- `Snakefile.py`: Core pipeline logic with Snakemake rules (simulation, analysis, docking prep, etc.)

## Dependencies

Snakemake

GROMACS

Python 3.8+

NumPy, RDKit, MDTraj, etc.


## Add a Conda Environment (Optional but Helpful)

You can create an environment.yaml:

name: md_pipeline
channels:
  - conda-forge
  - bioconda
dependencies:
  - snakemake
  - python=3.9
  - numpy
  - mdanalysis
  - rdkit
  - openmm

## 📊 Outputs

Clustered representative snapshots (.pdb)
MD trajectories and logs
Docking scores and poses (.sdf, .pdbqt)
Summary reports (.csv, visual plots of clustering/docking)


