# [CRISPR Analysis Software](https://github.com/ryandward/crispr_experiment)
A toolkit for analyzing CRISPR experiments

## Setup
1. Install dependencies:
   - Install [Miniforge](https://github.com/conda-forge/miniforge#Install) -- copy and paste the curl/wget script provided on the install page
2. Clone this repository:
   ```bash
   git clone https://github.com/ryandward/crispr_experiment.git
   cd crispr_experiment
   ```
3. Create and activate the environment:
   ```bash
   mamba env create -f environment.yml -n crispr_experiment
   mamba activate crispr_experiment
   ```

## Usage
Activate the environment before each session:
```bash
mamba activate crispr_experiment
```

## Updating
Pull the latest changes:
```bash
git pull
```

## GUIs
- extensible_GUI:
  - Main application entry point that provides access to all modules.
  - Features "Relaunch Application" and "Update Software" options.
  - Access all tools through a unified, user-friendly interface.

Main Tools:
- find_guides_gui:
  - Design new CRISPR guides from any genome sequence.
  - Select parameters like barcode length and PAM sequence.
  - Generates guides based on your specifications.

- design_mismatches_gui:
  - Create barcode variants by inserting mismatches into guides.
  - Check for potential off-targets in your genome.
  - Filter out problematic or overly similar guides.

Extra Tools:
- assembly_finder_gui:
  - Search and download genome assemblies from NCBI.
  - Input accession numbers to fetch sequences.
  - Track and manage downloaded assemblies.

[Guide](https://github.com/ryandward/crispr_experiment) | 
[Issues](https://github.com/ryandward/crispr_experiment/issues) | 
[Contributing](https://github.com/ryandward/crispr_experiment/blob/main/CONTRIBUTING.md)
