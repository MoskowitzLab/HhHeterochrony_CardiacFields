import pandas as pd
import scanpy as sc
import anndata as ad


def format_adata(
    adata_dir: str,
    adata_fname: str,
    metadata_path: str
) -> ad.AnnData:
    """Loads Pijuan dataset Anndata, adds metadata for
    plotting (original UMAP coordinates and cell type colors
    from paper), removes "mixed gastrulation" cells, and
    normalizes data for plotting. Saves and returns adata 
    object.
    
    Parameters
    ----------
    adata_dir
        a string defining the path to the directory 
        containing the Pijuan anndata object
        ***needs to end in /***

    adata_fname
        a string defining file name of the Pijuan anndata
        object 
    
    metadata_path
        a string defining the path to the Pijuan metadata
        file (contains UMAP colors and UMAP coordinates)

    Returns
    -------
    adata
        an anndata object containing the formatted Pijuan
        data"""

    adata = sc.read_h5ad(adata_dir+adata_fname)
    metadata = pd.read_csv(metadata_path)
    metadata.dropna(inplace=True)

    umap_coordinates = metadata[["UMAP_1", "UMAP_2"]].to_numpy()
    adata.obsm["UMAP_coordinates"] = umap_coordinates

    # color metadata is per cell-- map color for each cell
    # to its corresponding cell type and remove duplicates
    # to get list of unique cell types and corresponding
    # colors
    metadata["color"] = "#" + metadata["colour"].astype(str)
    colors = metadata["color"].copy().reset_index(drop=True)
    celltypes = adata.obs[
        "celltype_PijuanSala2019"].copy().reset_index(drop=True)
    color_mapping_df = pd.concat([celltypes, colors], axis=1)
    color_mapping_df.drop_duplicates(inplace=True)
    color_dict = dict(zip(
        color_mapping_df["celltype_PijuanSala2019"],
        color_mapping_df["color"]))
    
    adata.uns["color_dict"] = color_dict

    mask = adata.obs["stage"] != "Mixed gastrulation"
    adata = adata[mask, :]

    adata.layers["raw_counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.layers["norm_counts"] = adata.X.copy()
    adata.raw = adata

    adata.write_h5ad(adata_dir+"Pijuan_adata_formatted.h5ad")

    return adata


def filter_adata(
    adata: ad.AnnData, 
    adata_dir: str
) -> ad.AnnData:
    """Removes cell type cluters for stages at
    which there are fewer than 20 cells. Saves 
    and returns adata object.
    
    Parameters
    ----------
    adata
        an anndata object containing the formatted Pijuan
        data

    adata_dir
        a string defining the path to the directory 
        containing the Pijuan anndata object
        ***needs to end in /***

    Returns
    -------
    adata_filtered
        an anndata object containing the filtered Pijuan
        data"""

    #get cell types and stages from adata object
    celltypes = adata.obs["celltype_PijuanSala2019"]
    stages = adata.obs["stage"]

    #create df from adata X with celltype and stages columns
    X = adata.raw.X.toarray()
    df_ps2019 = pd.DataFrame(X, columns=adata.raw.var_names)
    df_ps2019["celltype"] = celltypes.values

    #make stages ordered categorical to maintain order for later plotting
    df_ps2019["stage"] = pd.Categorical(
        stages.values,
        categories=[
            "E6.5", "E6.75", "E7.0", "E7.25", "E7.5", "E7.75", "E8.0", "E8.25", "E8.5"
        ],
        ordered=True
    )

    #group df by celltype and stage to compute the total 
    #number of cells for each cell type for each stage
    #(creates a summary df with this info)
    df_celltypes = df_ps2019.groupby(
        ["celltype", "stage"], observed=True).agg(
            n_cells=("celltype", "count")).reset_index()

    #keep only the celltypes/stages combos where n_cells
    #is at least 20
    df_threshold = df_celltypes.copy(deep=True)[
        df_celltypes["n_cells"] >= 20]

    #format the celltype/stage combos as tuples to set 
    #up mask for adata
    keep = df_threshold[["celltype", "stage"]]
    keep_tuples = [tuple(x) for x in keep.values]

    #mask rows in adata that match both celltype and stage
    mask = adata.obs.apply(lambda row:
                           (row["celltype_PijuanSala2019"],
                            row["stage"]) in keep_tuples, 
                            axis=1)
    adata_filtered = adata[mask, :].copy()

    #save filtered adata
    adata_filtered.write_h5ad(adata_dir+"Pijuan_adata_filtered.h5ad")

    return adata_filtered