# MD Simulation & Ensemble Docking Workflow

This repository contains a modular and automated Snakemake pipeline for conducting **Molecular Dynamics (MD) simulations**, followed by **ensemble docking** using GROMACS, OpenMM, and docking tools such as Glide or AutoDock Vina. This pipeline was developed to explore protein–ligand interactions across dynamic conformations and is especially suited for high-throughput virtual screening efforts.

## 📌 Features

- Prepares protein structures for MD simulations.
- Performs energy minimization, equilibration, and production runs using GROMACS.
- Aligns and extracts representative MD snapshots using PCA and KMeans clustering.
- Prepares multiple receptor conformations for ensemble docking.
- Executes automated docking using Vina/Glide.
- Compatible with GPU-enabled HPC environments using SLURM.

## ⚙️ Requirements

- [Snakemake ≥ 7.0](https://snakemake.readthedocs.io)
- Conda (with Mamba recommended)
- GROMACS
- OpenMM
- Python ≥ 3.8 (with NumPy, MDAnalysis, MDTraj, scikit-learn, RDKit, etc.)
- AutoDock Vina or Glide (license required for Glide)
- SLURM (for HPC scheduling)

## Method Summary

**MD Simulations**: Conducted using GROMACS or OpenMM with system preparation (solvation, neutralization), followed by NPT equilibration and a production run.
Snapshot Extraction: Principal Component Analysis (PCA) followed by KMeans clustering is used to extract diverse representative structures.
Docking: Docking is performed against all representative receptor conformations to account for dynamic flexibility ("ensemble docking").

## 📊 Outputs

Clustered representative snapshots (.pdb)
MD trajectories and logs
Docking scores and poses (.sdf, .pdbqt)
Summary reports (.csv, visual plots of clustering/docking)


