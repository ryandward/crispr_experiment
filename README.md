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
   mamba create --name crispr_experiment --file spec-file.txt
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
mamba create --name crispr_experiment --file spec-file.txt
```

## Development
After adding new packages, update the environment files:
```bash
mamba list --explicit | grep -v "#" | grep -v "file:" > spec-file.txt
mamba env export | grep -v "prefix" > environment.yml

git add environment.yml spec-file.txt
git commit -m "Updated environment with new packages"
git push
```

[Guide](https://github.com/ryandward/crispr_experiment) | 
[Issues](https://github.com/ryandward/crispr_experiment/issues) | 
[Contributing](https://github.com/ryandward/crispr_experiment/blob/main/CONTRIBUTING.md)
