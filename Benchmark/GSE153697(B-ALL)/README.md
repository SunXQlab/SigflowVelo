# SigflowVelo is benchmarked against other methods using the relapsed B-cell acute lymphoblastic leukemia (B-ALL) dataset GSE153697

# Introduction
Here we introduce the detailed process of re-implementation of SigflowVelo against other 
methods including TFVelo, scVelo, UniTVelo and DeepVelo.  

The project follows a standardized directory structure:
- `data/`: Contains all necessary data files
- `environment/`: Stores environment configuration files
- `scripts/`: Houses all code files
- `results/`: Output directory for saving execution results

## Workflow

### 1. Clone the repository from GitHub
```bash
git clone https://github.com/SunXQlab/SigflowVelo.git
```

### 2. Reconstruct the environment
To ensure exact reproducibility, please recreate the virtual environment using the provided 
`method_environment.yml` file and replace `method` with the exact name of the method.  

Take SigflowVelo as an example.

```bash
cd SigflowVelo/benchmark/GSE153697(B-ALL)/environment
conda env create -f SigflowVelo_environment.yml
```

### 3. Running SigflowVelo
All the scripts of SigflowVelo are stored in `SigflowVelo/benchmark/GSE153697(B-ALL)/scripts/SigflowVelo/`
- `1_fsvelo_data_prepare.py`: The preparation of SigflowVelo's input.
- `2_Malignant_cells_imputate.R`:Imputation of missing values in datasets.
- `3_SigflowVelo.py`:Model training and velocity inference of SigflowVelo.
- `4_plot.py`:Visualizations of the velocity streamline and latent time inferred by SigflowVelo, corresponding to Fig.3 g.
- `5_validity.py`:Metric calculations for comparison with other methods and various visualizations, corresponding to Fig.3 e.

### 4. Running TFVelo
All the scripts of TFVelo are stored in `SigflowVelo/benchmark/GSE153697(B-ALL)/scripts/TFVelo/`. To install TFVelo, please clone 
the GitHub repository to your local machine: `https://github.com/xiaoyeye/TFvelo.git`. Then
install TFVelo according to its instructions.
- `1_TFvelo_run_demo.py`: Model training of TFVelo.
- `2_TFvelo_analysis_demo.py`: Velocity inference of TFVelo.
- `3_plot.py`: Various visualizations of TFVelo, corresponding to Fig.3 e and Supplement Fig.5 a-b.

### 5. Running scVelo
All the scripts of SigflowVelo are stored in `SigflowVelo/benchmark/GSE153697(B-ALL)/scripts/scVelo/`  
- `1_scVelo.py`: Velocity inference of scVelo.
- `2_plot.py`: Various visualizations of scVelo, corresponding to Fig.3 e and Supplement Fig.5 a-b.

### 6.Running UniTVelo
All the scripts of SigflowVelo are stored in `SigflowVelo/benchmark/GSE153697(B-ALL)/scripts/UniTVelo/`  
- `1_UniTVelo.py`: Velocity inference of UniTVelo.
- `2_plot.py`: Various visualizations of UniTVelo, corresponding to Fig.3 e and Supplement Fig.5 a-b.

### 7.Running DeepVelo
All the scripts of SigflowVelo are stored in `SigflowVelo/benchmark/GSE153697(B-ALL)/scripts/DeepVelo/`  
- `1_DeepVelo.py`: Velocity inference of DeepVelo.
- `2_plot.py`: Various visualizations of DeepVelo, corresponding to Fig.3 e and Supplement Fig.5 a-b.