# SigflowVelo is benchmarked against other methods on the multiple myeloma (MM) dataset GSE164551

# Introduction
Here we introduce the detailed process of re-implementation of SigflowVelo against TFVelo. 

The project follows a standardized directory structure:
- `data/`: Contains all necessary data files
- `environment/`: Stores environment configuration files
- `scripts/`: Houses all code files

## Workflow

### 1. Clone the repository from GitHub
```bash
git clone https://github.com/SunXQlab/SigflowVelo.git
```

### 2. Reconstruct the environment
To ensure exact reproducibility, please create virtual environments using Python's built-in venv module. Different methods require different Python versions: SigflowVelo requires Python 3.10.8, and TFVelo requires Python 3.8.10. The dependencies for each method are listed in the corresponding method_requirements.txt files (e.g., SigflowVelo_requirements.txt, TFVelo_requirements.txt). 

### 3. Running SigflowVelo
All the scripts of SigflowVelo are stored in `SigflowVelo/benchmark/GSE164551(MM)/scripts/SigflowVelo/`
- `1_fsvelo_data_prepare.py`: The preparation of SigflowVelo's input.
- `2_process_realdata.py`: Perform the first stage of preprocessing on the dataset
- `3_imputate.R`:Imputation of missing values in datasets.
- `4_SigflowVelo.py`:Model training and velocity inference using SigflowVelo.
- `5_plot_latent_time_stream.py`:Visualizations of the velocity streamline and latent time inferred by SigflowVelo.
- `6_Cos_Similarity.py`:Evaluates the directionality of inferred velocity by calculating cosine similarity and visualizing the results as a rose plot.
- `7_KS_statistic.py`:Assesses the separation of latent time distributions between conditions using the Kolmogorov–Smirnov statistic.

### 4. Running TFVelo
All the scripts of TFVelo are stored in `SigflowVelo/benchmark/GSE164551(MM)/scripts/TFVelo/`. To install TFVelo, please clone 
the GitHub repository to your local machine: `https://github.com/xiaoyeye/TFvelo.git`. Then install TFVelo according to its instructions.
- `1_data_processed.R`: Preprocessing of the dataset
- `2_TFvelo_run.py`: Model training of TFVelo.
- `3_TFvelo_analysis.py`: Velocity inference of TFVelo.
- `4_merged_scvelo_three_svgs.py`: Various visualizations of TFVelo.
- `5_Cos_Similarity.py`:Evaluates the directionality of inferred velocity by calculating cosine similarity and visualizing the results as a rose plot.
- `6_KS_statistic.py`:Assesses the separation of latent time distributions between conditions using the Kolmogorov–Smirnov statistic.