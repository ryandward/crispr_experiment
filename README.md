# CRISPR Analysis Software Setup
## First Time Setup
1. Install required software:
   - Install Miniforge from https://github.com/conda-forge/miniforge
2. Get the analysis tools:
   ```
   git clone https://github.com/ryandward/crispr_experiment.git
   cd crispr_experiment
   ```
3. Set up your software environment:
   ```
   mamba create --name crispr_experiment --file spec-file.txt
   mamba activate crispr_experiment
   ```
## Daily Use
Every time you start working, activate your environment:
```
mamba activate crispr_experiment
```
## Getting Updates
When new versions are released:
```
git pull
mamba create --name crispr_experiment --file spec-file.txt
```
## Creating Updates
After installing new packages, save both strict and flexible versions:
```
mamba list --explicit | grep -v "#" | grep -v "file:" > spec-file.txt
mamba env export | grep -v "prefix" > environment.yml
```
```
git add environment.yml spec-file.txt
git commit -m "Updated environment with new packages"
git push
```
