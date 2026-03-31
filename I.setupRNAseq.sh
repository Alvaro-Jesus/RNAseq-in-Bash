#!/usr/bin/env bash
set -euo pipefail

echo "    This is a script to perform automatically the"
echo "    INSTALLING MINICONDA + RNA-SEQ ENV by Alvaro Huamani Ortiz        "
echo "    Version 1.0, December 2nd 2025, 6 steps ~~discoverthevibes~~"

# -----------------------------
# 1. DOWNLOAD AND INSTALL MINICONDA
# -----------------------------
echo "STEP 1.  Downloading Miniconda installer..."

# Choose installer depending on OS:
INSTALLER="Miniconda3-latest-Linux-x86_64.sh"

curl -LO https://repo.anaconda.com/miniconda/$INSTALLER

echo "STEP 1 DONE. Miniconda download complete"

echo "STEP 2. Installing Miniconda silently, no clicks!!!"
bash $INSTALLER -b -p "$HOME/miniconda3"

# Initialize conda
echo "Initializing conda in the same shell..."
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"

# Evitar que base se active automáticamente
conda config --set auto_activate_base false

echo "STEP 2 DONE. Conda initialized, all good"
# -----------------------------
# 2. UPDATE CONDA + INSTALL MAMBA
# -----------------------------
echo "STEP 3.  Updating conda and installing mamba..."

conda update -n base -y conda
conda install -n base -c conda-forge -y mamba

echo  "STEP 3 is done."
# -----------------------------
# 3. CONFIGURE CHANNELS
# -----------------------------
echo "STEP 4. Configuring conda channels..."

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict


echo "STEP 4 is done. Almost there, dude!"
# -----------------------------
# 4. CREATE RNA-SEQ ENVIRONMENT
# -----------------------------
echo "STEP 5.  Creating environment: rnaseq ..."

mamba create -n rnaseq -y \
    sra-tools \
    hisat2 \
    fastp \
    multiqc \
    samtools \
    fastqc \
    pigz \
    seqtk \
    parallel

echo "STEP 5 is done. One more step and is over"
# -----------------------------
# 5. VERIFY INSTALLATION
# -----------------------------
echo "STEP 6. Verifying tools inside environment..."

conda activate rnaseq
prefetch --version || true
fasterq-dump --version || true
hisat2 --version
fastp --version
multiqc --version
samtools --version
fastqc --version
seqtk 2>&1 | head -n 1 || true
parallel --version | head -n 1

echo "DONE! Conda Environment 'rnaseq' is ready. Time to try the next script!"
