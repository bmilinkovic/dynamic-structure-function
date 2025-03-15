# Phase Coherence Analysis Scripts

This directory contains optimized scripts for analyzing phase coherence in fMRI data across different drug conditions and consciousness states. The analysis is implemented in both MATLAB and Python.

Author: Borjan Milinkovic

Updated: 2025

## MATLAB Implementation

### Scripts Overview

1. `kmeans_analysis_hilbert_example_KETA_psy.m`

   - Analyzes phase coherence in Ketamine vs Placebo conditions

2. `kmeans_analysis_hilbert_example_LSD_enzo.m`

   - Analyzes phase coherence in LSD vs Placebo conditions

3. `kmeans_analysis_hilbert_example_MDMA_enzo.m`
   - Analyzes phase coherence in MDMA vs Placebo conditions

### MATLAB Requirements

- MATLAB (tested on R2019b or later)
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

## Python Implementation

### Jupyter Notebook

`phase_coherence_analysis_documented.ipynb` provides a Python implementation of the phase coherence analysis, with the following features:

1. **Data Loading and Preprocessing**

   - Loading BOLD time series data
   - Z-score normalization
   - Signal preprocessing

2. **Phase Analysis**

   - Hilbert transform for phase extraction
   - Phase synchronization computation
   - Pattern extraction

3. **K-means Clustering**

   - Identification of brain states
   - Cluster analysis
   - State probability computation

4. **Visualization**
   - Brain state patterns
   - State probabilities
   - Temporal dynamics
   - Transition matrices

### Python Requirements

See `requirements.txt` for specific versions. Main dependencies:

- Python 3.8+
- NumPy
- Pandas
- SciPy
- Matplotlib
- Seaborn
- Scikit-learn
- Networkx

## Input Data Requirements

Each analysis requires:

1. Drug-specific time series data (e.g., `KETA_TS_FC.mat`, `LSD_TS_FC.mat`, or `MDMA_TS_FC.mat`)
2. Structural connectivity data (`Structural.mat`)

## Usage

### MATLAB Scripts

1. Ensure required data files are in your MATLAB path
2. Set your working directory to the script location
3. Run the desired script

### Python Notebook

1. Install requirements:
   ```bash
   pip install -r requirements.txt
   ```
2. Launch Jupyter:
   ```bash
   jupyter notebook
   ```
3. Open `phase_coherence_analysis_documented.ipynb`
4. Follow the step-by-step analysis with detailed explanations

## Output

Both implementations generate:

1. Brain state patterns identified through k-means clustering
2. Correlation between states and structural connectivity
3. State probability distributions
4. Temporal dynamics visualizations

## Notes

- The scripts use parallel processing for k-means clustering when available
- Outlier removal is performed using cityblock distance
- Phase patterns are computed using the Hilbert transform
- Data is preprocessed using z-scoring and detrending

## Optimizations Made

1. **Improved Code Organization**

   - Clear section headers with detailed comments
   - Consistent variable naming conventions
   - Modular code structure
   - Proper error handling

2. **Enhanced Documentation**

   - Detailed headers explaining purpose and requirements
   - Inline comments explaining complex operations
   - Clear parameter descriptions
   - Documentation of expected inputs and outputs

3. **Code Efficiency**

   - Pre-allocation of arrays
   - Vectorized operations where possible
   - Optimized data concatenation
   - Improved memory management

4. **Visualization Improvements**
   - Consistent figure layouts
   - Better subplot organization
   - Added colorbars
   - Improved titles and labels
