# CRISPR Experiment Environment Setup

## Repository Setup

To clone this repository:
```
git clone https://github.com/ryandward/crispr_experiment.git
```
To pull updates:
```
git pull
```

## Environment Management

### Create the environment
```
mamba create --name crispr_experiment --file spec-file.txt
conda activate crispr_experiment
```

### Export or update
```
mamba env export | grep -v "prefix" > environment.yml
mamba env update --file environment.yml --prune
```

### Strict version control
```
mamba list --explicit | grep -v "#" | grep -v "file:" > spec-file.txt
```

## Detailed Setup for Beginners

1. Install Miniconda or Anaconda from the official website (https://docs.conda.io/en/latest/miniconda.html).  
2. Once installed, open a terminal (or Anaconda Prompt on Windows) and install Mamba:
   ```
   conda install mamba -n base -c conda-forge
   ```
3. Clone the repository:
   ```
   git clone https://github.com/ryandward/crispr_experiment.git
   ```
   Or, if you prefer not to use Git, download and unpack the ZIP from GitHub.

## Setting Up the Environment

1. Navigate to the cloned folder:
   ```
   cd crispr_experiment
   ```
2. Create a new environment using the spec file:
   ```
   mamba create --name crispr_experiment --file spec-file.txt
   ```
3. Activate the environment:
   ```
   conda activate crispr_experiment
   ```

## Updating the Environment

- To get the latest code changes:
  ```
  git pull
  ```
- To update the existing environment after changes:
  ```
  mamba env update --file environment.yml --prune
  ```# crispr_experiment
