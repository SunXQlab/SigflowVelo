# SigflowVelo: A physics-informed framework for reconstructing signaling-driven cell state transitions and drug resistance from single-cell data

![](./SigflowVelo.png)

SigflowVelo is a physics-informed computational framework designed to reconstruct cell-state transition dynamics driven by intercellular signaling flow. By integrating single-cell RNA sequencing (scRNA-seq) snapshots with a kinetic model, it couples cell–cell communication (CCC) directly with continuous phenotypic plasticity.

## Installation

SigflowVelo can be installed directly via `pip`.

```python
conda create -n sigflowvelo python=3.11
conda activate sigflowvelo
pip install SigflowVelo==0.1.1
```

## Data Preparation

Within the `SigflowVelo_analysis` folder, `1_PDAC_CCC.R`, `2_flowsig.py`, and `3_fsvelo_data_prepare.py` are primarily used for single-cell data preprocessing, the construction of cell-cell communication and regulatory networks, and formatting the model input matrix (`AnnData`). The `4_process_realdata.py` script utilizes the prepared data for neural network training and the final inference of transcriptional dynamics velocity.
