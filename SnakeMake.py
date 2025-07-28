import os
from pathlib import Path
import glob

configfile:
    "config.yaml"

PATH = "/home/naqme/md-pipeline/"
SCRIPTS = PATH + "scripts"

# Define parameters for .mdp (Molecular Dynamics Parameter) files from the config file
# The .mdp files are stored in md-pipeline/gromacs-utils

NVT_TIME = config["nvt_time"]
NVT_DT = config["nvt_dt"]

NPT_TIME = config["npt_time"]
NPT_DT = config["npt_dt"]

MD_TIME = config["md_time"]
MD_DT = config["md_dt"]
SAVE_FREQ = config["save_freq"]

# Configuration parameters from the config file

PDB_FILE = config["pdb_file"]
PROTEIN_NAME = config["protein_name"]
OUTPUT_PATH = config["output_path"]
UTILS_FILES = config["utils_files"]
BOX_SIZE = config["box_size"]
BOX_SHAPE = config["box_shape"]
ION_CONCENTRATION = config["ion_concentration"]
MD_LENGTH = config["simulation_length"]

# PCA and KMeans Analysis configuration

SELECTION = config["selection"]
COMPONENTS = config["pca_components"]
ML = PATH + "ML_files"
XVG = PATH + "xvg_files"
os.makedirs(XVG, exist_ok=True)


rule all:
    input:
        "end.txt"


# ===================
# Step initializing: remove the simulation parameter files (.mdp files) for a new simulation
# ===================
# Remove the parameter files.
rule init_cleanup:
    output:
        ".done"
    run:        
        param_files = [
            os.path.join(UTILS_FILES, "nvt.mdp"),
            os.path.join(UTILS_FILES, "npt.mdp"),
            os.path.join(UTILS_FILES, "md.mdp"),
            os.path.join(OUTPUT_PATH, "truncation_point.txt")
        ]

        for param_file in param_files:
            if os.path.isfile(param_file):
                os.remove(param_file)
                print(f"The existing {param_file} parameter file was removed.")
            else:
                print(f"No {param_file} was found! Moving on!")
                pass
        shell("touch .done")


# ===================
# Step -1: Generate/Refine the simulation parameter files (.mdp files) 
# ===================
# This rule generates the .mdp files for the NVT, NPT, and MD stages.
# It ensures that key simulation parameters (e.g., duration, time step, and save frequency) 
# are customizable based on user-defined configurations.
rule generate_parameter_files:
    input:
        ".done"
    output:
        UTILS_FILES + "/nvt.mdp",
        UTILS_FILES + "/npt.mdp",
        UTILS_FILES + "/md.mdp"
    params:
        nvtt = NVT_TIME,
        nvtd = NVT_DT,
        nptt = NPT_TIME,
        nptd = NPT_DT,
        mdt = MD_TIME,
        mdd = MD_DT,
        mdfreq = SAVE_FREQ
    shell:
        """
        python {SCRIPTS}/setting.py nvt-setting --duration {params.nvtt} --dt {params.nvtd} --o {output[0]} &&
        python {SCRIPTS}/setting.py npt-setting --duration {params.nptt} --dt {params.nptd} --o {output[1]} &&
        python {SCRIPTS}/setting.py md-setting --duration {params.mdt} --dt {params.mdd} --save_freq {params.mdfreq} --o {output[2]}
        """


# ===================
# Step 0: Prepare Protein
# ===================
# Add missing residues/atoms and hydrogen at pH 7.0
rule preprocessing:
    input:
        PDB_FILE
    output:
        OUTPUT_PATH + "/preprocessed" + PROTEIN_NAME + ".pdb"
    conda: 
        "envs.yaml"
    shell:
        'python {SCRIPTS}/pdbfixer_script.py {input[0]} {output[0]}' 


# ===================
# Step 1: Generate Topology
# ===================
# Generate a GROMACS topology file and initial GRO file from the input PDB file
rule generate_topology:
    input:
        OUTPUT_PATH + "/preprocessed" + PROTEIN_NAME + ".pdb",
        UTILS_FILES + "/posre.itp"
    output:
        OUTPUT_PATH + "/file_processed_" + PROTEIN_NAME + ".gro",
        OUTPUT_PATH + "/topol_" + PROTEIN_NAME + ".top"
    singularity: '../singularity/gromacs.sif'
    shell:
        """
        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\
        gmx pdb2gmx -f {input[0]} -o {output[0]} -water tip3p -p {output[1]} -i {input[1]} -ff charmm36-jul2022 -ter -ignh << EOF
        0
        0
        EOF"
        """

# ============================
# Step 2: Define Simulation Box for solvation and ionization steps
# ============================
# Generate the simulation box, including the water molecules and ions (to neutralise the system and add ionic strength)
rule solvation:
    input:
        OUTPUT_PATH + "/file_processed_" + PROTEIN_NAME + ".gro",
        OUTPUT_PATH + "/topol_" + PROTEIN_NAME + ".top"
    params:
        box_size = BOX_SIZE,
        box_shape= BOX_SHAPE
    output:
        OUTPUT_PATH + "/file_newbox_" + PROTEIN_NAME + ".gro",
        OUTPUT_PATH + "/file_solv_" + PROTEIN_NAME + ".gro"
    singularity: '../singularity/gromacs.sif'
    shell:
        """
        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\
        
        # Step 1
        gmx editconf -f {input[0]} -o {output[0]} -c -d {params.box_size} -bt {params.box_shape} &&\

        # Step 2
        gmx solvate -cp {output[0]} -cs spc216.gro -o {output[1]} -p {input[1]}"

        """

rule ionization:
    input:
        UTILS_FILES + "/ions.mdp",
        OUTPUT_PATH + "/file_solv_" + PROTEIN_NAME + ".gro",
        OUTPUT_PATH + "/topol_" + PROTEIN_NAME + ".top"
    output:
        OUTPUT_PATH + "/ions_" + PROTEIN_NAME + ".tpr",
        OUTPUT_PATH + "/file_solv_ions_" + PROTEIN_NAME + ".gro"
    params:
        conc = ION_CONCENTRATION
    singularity: '../singularity/gromacs.sif'
    shell:
        """
        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\
        # Step 1
        gmx grompp -f {input[0]} -c {input[1]} -p {input[2]} -o {output[0]} &&\

        # Step 2
        gmx genion -s {output[0]} -o {output[1]} -p {input[2]} -pname NA -nname CL -neutral -conc {params.conc} << EOF
        13
        EOF"
        """

# ============================
# Step 3: Energy Minimization
# ============================
# Process and run energy minimization
rule energy_minimization:
    input:
        UTILS_FILES + "/minim.mdp",
        OUTPUT_PATH + "/file_solv_ions_" + PROTEIN_NAME + ".gro",
        OUTPUT_PATH + "/topol_" + PROTEIN_NAME + ".top"
    output:
        OUTPUT_PATH + "/em_" + PROTEIN_NAME  + ".tpr",
        OUTPUT_PATH + "/em_" + PROTEIN_NAME  + ".gro",
        XVG + "/potential_" + PROTEIN_NAME + ".xvg"
    singularity: '../singularity/gromacs.sif'
    shell:
        """
        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\
        # Step 1
        gmx grompp -f {input[0]} -c {input[1]} -p {input[2]} -o {output[0]} &&\

        # Step 2
        gmx mdrun -v -s {output[0]} -c {output[1]} &&\

        # Step 3
        gmx energy -s {output[0]} -o {output[2]} << EOF
        10
        0
        EOF"
        """

# ============================
# Step 4: NVT Equilibration
# ============================
# Process and run NVT (constant volume and temperature) equilibration
rule equilibration_NVT:
    input:
        UTILS_FILES + "/nvt.mdp",
        OUTPUT_PATH + "/em_" + PROTEIN_NAME + ".gro",
        OUTPUT_PATH + "/topol_" + PROTEIN_NAME + ".top"
    output:
        OUTPUT_PATH + "/nvt_" + PROTEIN_NAME + ".tpr",
        OUTPUT_PATH + "/nvt_" + PROTEIN_NAME + ".gro",
        OUTPUT_PATH + "/nvt_" + PROTEIN_NAME + ".cpt",
        OUTPUT_PATH + "/nvt_" + PROTEIN_NAME + ".edr",
        XVG + "/temprature_" + PROTEIN_NAME + ".xvg"
    singularity: '../singularity/gromacs.sif'
    shell:
        """
        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\
        # Step 1
        gmx grompp -f {input[0]} -c {input[1]} -r {input[1]} -p {input[2]} -o {output[0]} &&\

        # Step 2
        gmx mdrun -v -s {output[0]} -c {output[1]} -cpo {output[2]} -e {output[3]} &&\

        # Step 3
        gmx energy -f {output[3]} -o {output[4]} << EOF
        16
        0
        EOF"
        """


# ============================
# Step 5: NPT Equilibration
# ============================
# Process and run NPT (constant pressure and temperature) equilibration
rule equilibration_NPT:
    input:
        UTILS_FILES + "/npt.mdp",
        OUTPUT_PATH + "/nvt_" + PROTEIN_NAME + ".gro",
        OUTPUT_PATH + "/nvt_" + PROTEIN_NAME + ".cpt",
        OUTPUT_PATH + "/topol_" + PROTEIN_NAME + ".top"
    output:
        OUTPUT_PATH + "/npt_" + PROTEIN_NAME + ".tpr",
        OUTPUT_PATH + "/npt_" + PROTEIN_NAME + ".gro",
        OUTPUT_PATH + "/npt_" + PROTEIN_NAME + ".cpt",
        OUTPUT_PATH + "/npt_" + PROTEIN_NAME + ".edr",
        XVG + "/pressure_" + PROTEIN_NAME + ".xvg"
    singularity: '../singularity/gromacs.sif'
    shell:
        """
        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\
        # Step 1
        gmx grompp -f {input[0]} -c {input[1]} -r {input[1]} -p {input[3]} -o {output[0]} &&\

        # Step 2
        gmx mdrun -v -s {output[0]} -c {output[1]} -cpo {output[2]} -e {output[3]} &&\

        # Step 3
        gmx energy -f {output[3]} -o {output[4]} << EOF
        18
        0
        EOF"
        """
# ============================
# Step 6: MD Production
# ============================
# Run the main production run
rule production:
    input:
        UTILS_FILES + "/md.mdp",
        OUTPUT_PATH + "/npt_" + PROTEIN_NAME + ".gro",
        OUTPUT_PATH + "/npt_" + PROTEIN_NAME + ".cpt",
        OUTPUT_PATH + "/topol_" + PROTEIN_NAME + ".top"
    output:
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.gro",
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.cpt",
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.edr",
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.tpr",
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.log",
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.xtc"
    singularity: '../singularity/gromacs.sif'
    shell:
        """
        singularity exec --nv ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\
        # Step 1
        gmx grompp -f {input[0]} -c {input[1]} -t {input[2]} -p {input[3]} -o {output[3]} &&\
        # Step 2
        gmx mdrun -v -s {output[3]} -c {output[0]} -cpo {output[1]} -e {output[2]} -g {output[4]} -x {output[5]} -nb gpu"
        """


# ============================
# Step 7A: Initial Analysis
# ============================
# To account for trajectory periodicity, remove rorational and translational motions of the system, measure RMSD, RMSF and radius of gyration
rule correct_periodicity:
    input:
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.tpr",
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.xtc"
    output:
        OUTPUT_PATH + "/md_noPBC" + PROTEIN_NAME + MD_LENGTH + "ns.xtc"
    singularity: '../singularity/gromacs.sif'
    shell:
        """
        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\

        # Account for periodicity
        gmx trjconv -s {input[0]} -f {input[1]} -o {output[0]} -pbc mol -center << EOF
        1
        0
        EOF"
        """

# Remove rotational and translational motions
rule remove_rot_tran:
    input:
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.tpr",
        OUTPUT_PATH + "/md_noPBC" + PROTEIN_NAME + MD_LENGTH + "ns.xtc"
    output:
        OUTPUT_PATH + "/trajectory_" + PROTEIN_NAME + MD_LENGTH + "ns.xtc"
    singularity: '../singularity/gromacs.sif'
    shell:
        """
        # Remove center of mass motion
        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\

        gmx trjconv -fit rot+trans -s {input[0]} -f {input[1]} -o {output[0]} << EOF
        1
        0
        EOF"
        """

rule RMSD:
    input:
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.tpr",
        OUTPUT_PATH + "/trajectory_" + PROTEIN_NAME + MD_LENGTH + "ns.xtc"
    output:
        XVG + "/rmsd_" + PROTEIN_NAME + MD_LENGTH + "ns.xvg"
    singularity: '../singularity/gromacs.sif'
    shell:
        """
        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\
        gmx rms -s {input[0]} -f {input[1]} -o {output[0]} -tu ns << EOF
        4
        1
        EOF"
        """

# A text file is generated to set a point for trajectory cut (from the starting frame)
rule truncation_point:
    input:
        ".done"
    output:
        time= OUTPUT_PATH + "/truncation_point.txt"
    run:
        truncation_at = MD_TIME * 0.3 * 1000
        with open (output.time, "w") as writefile:
            writefile.write(f"{int(truncation_at)}\n")
            print(f"The text file for trajectory truncation at {truncation_at} (ps) is generated in {OUTPUT_PATH} path. Please note that the time is in picosecond.")


# This rule truncates the inital 30% of the trajectory (according to the text file generated in time_trunc rule) to access the more stablized conformations of the protein for upcoming analysis.
rule truncate_trajectory:
    input:
        ".done",
        tpr = OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.tpr",
        old_traj = OUTPUT_PATH + "/trajectory_" + PROTEIN_NAME + MD_LENGTH + "ns.xtc",
        truncation_time= OUTPUT_PATH + "/truncation_point.txt"
    output:
        new_traj= OUTPUT_PATH + "/trajectory_truncated" + PROTEIN_NAME + MD_LENGTH + "ns.xtc"
    singularity: '../singularity/gromacs.sif'
    shell:
       """
        time=$(cat {input.truncation_time}) &&

        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\

        gmx trjconv -s {input.tpr} -f {input.old_traj} -o {output.new_traj} -b $time << EOF
        0
        EOF"
        """


rule RMSF:
    input:
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.tpr",
        OUTPUT_PATH + "/trajectory_truncated" + PROTEIN_NAME + MD_LENGTH + "ns.xtc"
    output:
        XVG + "/rmsf_" + PROTEIN_NAME + MD_LENGTH + "ns.xvg"
    singularity: '../singularity/gromacs.sif'
    shell:
        """
        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\

        gmx rmsf -s {input[0]} -f {input[1]} -o {output[0]} -res << EOF
        4
        EOF"    
        """


rule radius_of_gyration:
    input:
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.tpr",
        OUTPUT_PATH + "/trajectory_truncated" + PROTEIN_NAME + MD_LENGTH + "ns.xtc"
    output:
        XVG + "/gyration_" + PROTEIN_NAME + MD_LENGTH + "ns.xvg"
    singularity:
        '../singularity/gromacs.sif'
    shell:
        """
        singularity exec ../singularity/gromacs.sif bash -c "source /usr/local/gromacs/bin/GMXRC &&\

        # Calculate radius of gyration
        gmx gyrate -s {input[0]} -f {input[1]} -o {output[0]} << EOF
        1
        EOF"
        """


# ============================
# Step 7B: Advanced Analysis
# ============================
# The PCA and KMeans are calculated and then the representative snapshot of each cluster is extracted
rule perform_PCA_KMeans:
    input:
        OUTPUT_PATH + "/trajectory_truncated" + PROTEIN_NAME + MD_LENGTH + "ns.xtc",
        OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.gro"
    output:
        ML + "/optimal_K.txt"
    params:
        selec = SELECTION,
        pca_com = COMPONENTS,
        output = ML
    shell:
        """
        mkdir -p {ML}
        python {SCRIPTS}/1-pca_and_kmeans.py --traj {input[0]} --coord {input[1]} --selection {params.selec} --pca_components {params.pca_com} --output_path {params.output}
        """
            

rule clustering_and_snapshot_extraction:
    input:
        traj = OUTPUT_PATH + "/trajectory_truncated" + PROTEIN_NAME + MD_LENGTH + "ns.xtc",
        coord = OUTPUT_PATH + "/md_" + PROTEIN_NAME + MD_LENGTH + "ns.gro",
        k_value = ML + "/optimal_K.txt"
    output:
        ML + "/clusters and centriods.png",
        ML + "/representative_snapshots.txt",
        glob.glob(ML + "/snapshots_*.txt")
    params:
        selec = SELECTION,
        pca_com = COMPONENTS,
        output = ML
    run:
        with open (input.k_value, "r") as f:
            k = f.read().strip()
            print(k)
        
        shell("""
        python {SCRIPTS}/2-clustering_and_rep_snapshots.py --traj {input.traj} --coord {input.coord} --selection {params.selec} --pca_components {params.pca_com} --clusters {k} --start 0 --output_path {params.output}
        
        """)

# ============================
# Step 8: .xvg Files Visualisation
# ============================
# The .xvg files from previous steps, namely (pressure, temperature, radius of gyration, rmsd and rmsf) are visualised. This ensures the accuracy of the simulation.
rule visualisation:
    input:
        gyr = XVG + "/gyration_" + PROTEIN_NAME + MD_LENGTH + "ns.xvg",
        pot = XVG + "/potential_" + PROTEIN_NAME + ".xvg",
        tem = XVG + "/temprature_" + PROTEIN_NAME + ".xvg",
        press = XVG + "/pressure_" + PROTEIN_NAME + ".xvg",
        rmsd= XVG + "/rmsd_" + PROTEIN_NAME + MD_LENGTH + "ns.xvg",
        rmsf= XVG + "/rmsf_" + PROTEIN_NAME + MD_LENGTH + "ns.xvg"
    output:
        gyr_plot = XVG + "/gyration_" + PROTEIN_NAME + MD_LENGTH + "ns.png",
        pot_plot = XVG + "/potential_" + PROTEIN_NAME + ".png",
        tem_plot =  XVG + "/temprature_" + PROTEIN_NAME + ".png",
        press_plot = XVG + "/pressure_" + PROTEIN_NAME + ".png",
        rmsd_plot = XVG + "/rmsd_" + PROTEIN_NAME + MD_LENGTH + "ns.png",
        rmsf_plot= XVG + "/rmsf_" + PROTEIN_NAME + MD_LENGTH + "ns.png"
    shell:
        """
        python {SCRIPTS}/visualisation.py rgvisual --input {input.gyr} --output {output.gyr_plot} &&
        python {SCRIPTS}/visualisation.py restvisual --potential {input.pot} --temperature {input.tem} --pressure {input.press} --rmsd {input.rmsd} --rmsf {input.rmsf} \
        --potential_out {output.pot_plot} --temperature_out {output.tem_plot} --pressure_out {output.press_plot} --rmsd_out {output.rmsd_plot} --rmsf_out {output.rmsf_plot}
        """

# ============================
# Step 9: Output generation for rule all
# ============================
# This steps only ensures the data analysis files (from the last steps) are successfully generated.
# The output of this rule is used by the 'all' rule for easier accessibility.
rule end:
    input:
        XVG + "/gyration_" + PROTEIN_NAME + MD_LENGTH + "ns.png",
        XVG + "/potential_" + PROTEIN_NAME + ".png",
        XVG + "/temprature_" + PROTEIN_NAME + ".png",
        XVG + "/pressure_" + PROTEIN_NAME + ".png",
        XVG + "/rmsd_" + PROTEIN_NAME + MD_LENGTH + "ns.png",
        XVG + "/rmsf_" + PROTEIN_NAME + MD_LENGTH + "ns.png",
        ML + "/clusters and centriods.png",
        ML + "/representative_snapshots.txt",
        glob.glob(ML + "/snapshots_*.txt")
    output:
        "end.txt"
    run:
        import time
        from datetime import date

        ending_date = date.today()
        ending_date_str= ending_date.strftime("%A %d. %B %Y")
        ending_date_str_time= ending_date.ctime()

        with open(output[0], 'w') as endwrite:
            endwrite.write(F"The simulation of the protein {PROTEIN_NAME} has on date and time: {ending_date_str_time} ended.")

# ============================
# Step 10: Clean-up task
# ============================
# The file generated after modifying the parameter files (.mdp) needs to be removed 
# to allow for future modifications of these files for subsequent simulations.
onsuccess:
    if os.path.isfile(".done"):
        os.remove(".done")
        print("Cleanup complete: .done removed.")