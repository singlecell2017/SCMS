# SCMS: Single Cell Mass Spectrometry

## Installation  	
#### install from PyPI

    python -m pip install --index-url https://test.pypi.org/simple/ --no-deps scms==0.0.4
    
#### install from GitHub

	git clone https://github.com/singlecell2017/SCMS.git
	cd SCMS
	python setup.py install
    
## Quick Start

### Command line

    scms.py -f INPUT_FILE -o OUTPUT_PATH
    
    INPUT_FILE: A file of single-cell mass spectrometry: each two columns are a group which contains 
    the nuclear to mass ratio (m/z) and its signal strength.
    OUTPUT_PATH: output path
    
    Note: Cell names should be like groupname-xx. For example: A549-1,A549-2, gefitinib-6,gefitinib-8..., 
    name before the first '-' will be considered as the group name, so the group names are A549 and 
    gefitinib. Otherwise, -g should be provided: -g group_name1 group_name2 group_name3...

#### other parameters 
* metafile: Pre-processed file for direct visualization. Default:None
* ppm_threshold_peak: Peak error threshold for peak selection in the same cell. Default:10
* ppm_threshold_cell: Peak error threshold for combining different cells. Default: 20
* decrease: Peak selection mode. Default:True
* peak: Whether select peak for each cell. Default:True
* filter: Peaks that appear in less than n% cells will be filtered. Default:0.5
* knn: Whether use knn for missing value imputation. Default:True
* n_neighbors: KNN algorithm parameter. Default:5
* method: Dimensional reduction methods, including PLS,PCA,UMAP,tSNE. Default:all
* q_value: Adjusted p_value threshhold for extracting differential expressed signal among different groups. Default:0.05
* log2fold: Log-foldchange threshhold for extracting differential expressed signal among different groups. Default:0.5
* seed: Random seed

#### Output
Output will be saved in the output folder including:
* **signal.csv**: Processed signal strength data after peak selection and coordinate alignment
* **m_z.csv**: Origin m/z and coordinate alignment m/z correspondence for each cell
* **signal_filtered.csv**: Filtered signal strength data, peaks appear in less than n% cells are filtered
* **m_z_filtered.csv**: Filtered m/z correspondence for each cell
* **signal_filtered_knn.csv**: Signal strength data after KNN imputation
* **PLS/tSNE/PCA/UMAP.csv**: Dimensionality reduction results by different methods
* **de.csv**: Differentially expressed m/z among different groups: mean: mean of signal strength, mean: mean of signal strength, median: median of signal strength, mode: mode of signal strength, group1: log-foldchange of signal strength of group1 to other groups, t/f_score: t_score or f_score, p_values: p_values, q_values: adjusted p_values
* **violinplot.png**: Violin plot correspond to signal.csv, and show the distribution of the number of m_z signals of all cells
* **PLS/tSNE/PCA/UMAP.pdf**: Dimensionality reduction embeddings by different methods
* **heatmap.png**: Heatmap of differentially expressed m_z among different groups

#### Help
Look for more usage of scms

	scms.py --help
