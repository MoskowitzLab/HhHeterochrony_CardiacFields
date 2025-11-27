# Hh ligand detection in the Pijuan Dataset
## Set up and installation
We recommend using a conda environment for dependency management. The environment.yml file can be used to create an environment with the necessary dependencies (specified versions indicate versions used to produce results reported in the paper). 

To create an environment with the .yml file, run the following command:

```bash
$ conda env create -f environment.yml
```

## Code overview
The run_Hh_ligand_analysis.ipynb jupyter notebook can be used to generate all of the figures shown in the manuscript, and imports the functions from the format_adata.py and plot_Hh_ligands.py.

format_adata.py contains functions to load the Pijuan anndata file, add necessary metadata for plotting, remove "mixed gastrulation" cells, and filter the data to remove cell type clusters for stages at which there are fewer than 20 cells. 
The format_adata() function is based on the specific anndata and metadata files we used here, and should be edited as needed to import the data and ensure that the desired UMAP coordinates and cell type colors are present. 

plot_Hh_ligands.py contains functions to create and save custom UMAP plots, Hh ligand feature plots, dot plots, and stacked violin plots, as shown in the manuscript. 


