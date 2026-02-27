# SigflowVelo is benchmarked against other methods using the  monocyte-macrophage population within Triple-Negative Breast Cancer (TNBC) tumors treated with PD-L1 blockade and paclitaxel

# Introduction
Here we introduce the detailed process of re-implementation of SigflowVelo against other 
methods including TFVelo, scVelo, UniTVelo and DeepVelo.  

The project follows a standardized directory structure:
- `data/`: Contains all necessary data files
- `scripts/`: Houses all code files
- `results/`: Output directory for saving execution results
- `visualization_and_evaluation/`: scripts for generating interpretaion of results and evaluations.

## Workflow

### 1. Clone the repository from GitHub
```bash
git clone https://github.com/SunXQlab/SigflowVelo.git
```

### 2. Reconstruct the environment
To ensure exact reproducibility, please recreate the virtual environment using the provided 
`environment.yml` file and replace `method` with the exact name of the method.  

Take SigflowVelo as an example.

```bash
cd SigflowVelo/benchmark/GSE169246(TNBC)/
conda env create -f environment.yml
```

### 3. Running SigflowVelo
All the scripts of SigflowVelo are stored in `SigflowVelo/benchmark/GSE153697(B-ALL)/scripts/SigflowVelo/`
- `1_flowsig.py`: data processing by Flowsig.
- `2_fsvelo_data_prepare.R`:The preparation of SigflowVelo's input.
- `3_impute.py`:Imputation of missing values in datasets.
- `4_SigflowVelo.py`:Model training and velocity inference of SigflowVelo.




### 5. Running scVelo
All the scripts of SigflowVelo are stored in `SigflowVelo/benchmark/GSE153697(B-ALL)/scripts/scVelo/`  
- `scVelo_bm.py`: Velocity inference of scVelo.


### 6.Running UniTVelo
All the scripts of SigflowVelo are stored in `SigflowVelo/benchmark/GSE153697(B-ALL)/scripts/UniTVelo/`  
- `UniTVelo_bm.py`: Velocity inference of UniTVelo.


### 7.Running DeepVelo
All the scripts of SigflowVelo are stored in `SigflowVelo/benchmark/GSE153697(B-ALL)/scripts/DeepVelo/`  
- `DeepVelo_bm.py`: Velocity inference of DeepVelo.
