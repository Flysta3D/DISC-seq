#!/bin/bash
# Install necessary dependencies for the DISC-Seq pipeline

echo "Installing system dependencies..."
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y python3-pip r-base r-base-dev

# CentOS/RHEL
# sudo yum install -y python3-pip R

echo "Installing Python packages..."
pip3 install numpy matplotlib regex

echo "Installing R packages..."
R -e "install.packages(c('Seurat', 'ggplot2', 'data.table', 'dplyr', 'reshape2', 'optparse', 'stringr', 'future'), repos='https://cloud.r-project.org')"

echo "Installing additional tools..."
# Install fastp
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
sudo mv fastp /usr/local/bin/

# Installing fastq-multx (From ea-utils)
# You should install ea-utils according to your system
# For example, on Ubuntu:
# Ubuntu: sudo apt-get install ea-utils
# Or download precompiled binaries and compile from source if necessary

echo "Installing umi_tools..."
pip3 install umi_tools

echo "All dependencies installed successfully!"