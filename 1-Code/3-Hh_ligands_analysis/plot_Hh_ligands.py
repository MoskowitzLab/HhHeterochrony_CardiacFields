import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import seaborn as sns
import anndata as ad
import matplotlib.colors as mplc

#needed to draw axes as arrows for UMAPs
from matplotlib.patches import FancyArrowPatch

#need to manually add Arial font- matplotlib can't find it in default location
from matplotlib.font_manager import fontManager
fontManager.addfont("./Arial.ttf")
plt.style.use("paper.mplstyle.py")



def plot_UMAP_custom(
        adata: ad.AnnData, 
        umap_coords_obsm: str, 
        palette: dict,
        plot_name: str,
        fig_size: tuple,
        fig_title: str,
        fig_output_path: str
):
    """Creates a custom UMAP plot using specified
    coordinates with axes represented as arrows.

    Parameters
    ----------
    adata
        an anndata object containing the formatted Pijuan
        data

    umap_coords_obsm
        a string defining the name of the .obsm array in the
        adata that contains the UMAP coordinates for plotting
        (these are added to obsm of the adata object in 
        format_adata() using imported metadata)
    
    palette
        a dict defining the color to use for each cell type

    plot_name
        a string defining the name of the plot for saving
        ***needs to end in .svg file extension***
    
    fig_size
        a tuple defining the figure size (width, height)

    fig_title
        a string defining the figure title (to appear on the
        plot)

    fig_output_path
        a string defining the desired output directory for
        the plot
        ***needs to end in /***"""
    
    celltypes = "celltype_PijuanSala2019"

    fig, ax = plt.subplots(figsize=fig_size)

    # Remove spines, ticks, and tick labels
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    sc.pl.embedding(adata, basis=umap_coords_obsm, color=celltypes,
                    palette=palette, size=1, title=fig_title, ax=ax, 
                    show=False, legend_fontsize="xx-small")
             
    #force layout/render before drawing arrows to get ax limits
    fig.canvas.draw()

    # # Get updated limits
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()

    #add filled triangle arrows
    arrowprops = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8",
                    color="black", linewidth=0.5, transform=ax.transData,
                    mutation_scale=3, clip_on=False)

    origin_X_offset = x0 - 0.01*(x1 - x0) #0.01 for (7,3) figsize all time points; 
    origin_Y_offset = y0 - 0.01*(y1 - y0) #0.0255 for (1.5,1.5) figsize individual time points

    arrow_x = FancyArrowPatch((origin_X_offset, y0), (x1, y0), **arrowprops)
    arrow_y = FancyArrowPatch((x0, origin_Y_offset), (x0, y1), **arrowprops)

    ax.add_patch(arrow_x)
    ax.add_patch(arrow_y)

    # Label UMAP axes just outside the final limits
    ax.text((x0 + x1)/2, origin_Y_offset*1.15, "UMAP 1", ha="center", va="center", fontsize=6)
    ax.text(origin_X_offset*1.1, (y0 + y1)/2, "UMAP 2", ha="center", va="center", fontsize=6, rotation=90)
    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.savefig(fig_output_path+plot_name, bbox_inches="tight")


def calculate_vmin_vmax(
    adata: ad.AnnData,
    Hh_source_list: list,
) -> tuple[float, float]:
    """Calculates the global minimum and maximum
    Hh ligand expression to use for the color bar
    in plot_Hh_source_feature_plot (to maintain
    consistent scaling across plots).

    Parameters
    ----------
    adata
        an anndata object containing the formatted Pijuan
        data
    Hh_source_list
        a list of strings defining the Hh ligands that will
        be plotted with the same scaling

    Returns
    -------
    vmin
        a float defining the global minimum of Hh ligand
        expression to use for the color bar (to maintain
        consistent scaling across plots)
    
    vmax
        a float defining the global maximum of Hh ligand
        expression to use for the color bar (to maintain
        consistent scaling across plots)"""
    
    # Compute global min/max for all three genes for consistent colormap scaling
    expr = adata[:, Hh_source_list].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    vmin = expr.min()
    vmax = expr.max()

    return vmin, vmax
    

def plot_Hh_source_feature_plot(
        adata: ad.AnnData,
        Hh_source: str,
        umap_coords_obsm: str,
        fig_size: tuple,
        vmin: float,
        vmax: float,
        fig_output_path: str
):
    """Creates a feature plot in UMAP space for the
    specified Hh source ligand, using the same format
    as plot_UMAP_custom(). Can be used to plot all stages
    combined, or separated by stage.

    Parameters
    ----------
    adata
        an anndata object containing the formatted Pijuan
        data

    Hh_source
        a string defining the Hh ligand to plot

    umap_coords_obsm
        a string defining the name of the .obsm array in the
        adata that contains the UMAP coordinates for plotting
        (these are added to obsm of the adata object in 
        format_adata() using imported metadata)
    
    fig_size
        a tuple defining the figure size (width, height)

    vmin
        a float defining the global minimum of Hh ligand
        expression to use for the color bar (to maintain
        consistent scaling across plots)
    
    vmax
        a float defining the global maximum of Hh ligand
        expression to use for the color bar (to maintain
        consistent scaling across plots)

    fig_output_path
        a string defining the desired output directory for
        the plot
        ***needs to end in /***"""

    plot_name = Hh_source+"_feature_plot.svg"
    fig, ax = plt.subplots(figsize=fig_size)

    # Remove spines, ticks, and tick labels
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    colors_list = ["#F0F921", "#FDC328", "#F89441", "#EB6A5C", "#D84B6A", "#B63679", "#8B2C80", "#5B2482", "#2B1585", "#0D0887" ]
    cmap = mplc.LinearSegmentedColormap.from_list("", colors_list)

    # if legend:
    sc.pl.embedding(adata, basis=umap_coords_obsm, color=Hh_source,
                    cmap=cmap, size=1, title=Hh_source +" expression", ax=ax,
                    show=False, legend_fontsize="xx-small", vmin=vmin, vmax=vmax)
    fig.canvas.draw()

    # # Get updated limits
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()

    #############
    #This block of code makes the axes arrows and adds axis labels manually
    #############
    #make arrow object
    arrowprops = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8",
                    color="black", linewidth=0.5, transform=ax.transData,
                    mutation_scale=3, clip_on=False)
    #add offset so arrows intersect at origin (manually tune)
    origin_X_offset = x0 - 0.0255*(x1 - x0) #0.01 for (7,3) figsize all time points; 
    origin_Y_offset = y0 - 0.0255*(y1 - y0) #0.0255 for (1.5,1.5) figsize individual time points
    #draw arrows
    arrow_x = FancyArrowPatch((origin_X_offset, y0), (x1, y0), **arrowprops)
    arrow_y = FancyArrowPatch((x0, origin_Y_offset), (x0, y1), **arrowprops)
    ax.add_patch(arrow_x)
    ax.add_patch(arrow_y)

    # Label UMAP axes just outside the final limits
    ax.text((x0 + x1)/2, origin_Y_offset*1.15, "UMAP 1", ha="center", va="center", fontsize=6)
    ax.text(origin_X_offset*1.1, (y0 + y1)/2, "UMAP 2", ha="center", va="center", fontsize=6, rotation=90)
    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.savefig(fig_output_path+plot_name, bbox_inches="tight")


def plot_Hh_source_dot_plot(
        adata: ad.AnnData,
        Hh_source_list: list,
        fig_height: float,
        fig_output_path: str,
        exp_threshold: float=None,
        annotate_frac_exp: bool=False
):
    """Creates a dot plot of Hh_source ligand
    expression across stages for each cell type,
    with optional expression threshold for plotting
    a given cell type. Makes a figure of subplots,
    one for each specified ligand.

    Parameters
    ----------
    adata
        an anndata object containing the filtered Pijuan
        data

    Hh_source_list
        a list of strings defining the Hh ligands to plot

    fig_height
        a float defining the height of the figure (the width
        is defined by the number of ligands/subplots in the
        figure)
    
    fig_output_path
        a string defining the desired output directory for
        the plot
        ***needs to end in /***

    exp_threshold
        a float defining the threshold of % cells in the cluster
        expressing the ligand at a given stage for plotting (as 
        a decimal). The cell type must have at least one time
        point meeting the threshold to be shown on the plot.
        Default is None (all cell types will be plotted).

    annotate_frac_exp
        a boolean defining whether to add annotations above
        each dot indicating the fraction of cells in the cluster
        expressing the ligand. Default is False."""

    celltypes = adata.obs["celltype_PijuanSala2019"]
    stages  = adata.obs["stage"]

    X = adata.raw.X.toarray()
    df = pd.DataFrame(X, columns=adata.raw.var_names)
    df["celltype"] = celltypes.values
    df["stage"] = pd.Categorical(
        stages.values,
        categories=[
            "E6.5", "E6.75", "E7.0", "E7.25", "E7.5", "E7.75", "E8.0", "E8.25", "E8.5"
        ],
        ordered=True
    )
    #############
    # This code block is necessary to get the limits for mean expression and % expression
    # for the plot.
    #############
    # Compute mean and percent expressing for each gene in Hh_source_list, then get global min/max
    # across lignads.
    mean_expr_list = []
    pct_expr_list = []
    if exp_threshold is not None:
        for gene in Hh_source_list:
            #df_gene consists of a row for each cell, and columns with cell type
            #annotation, stage, and gene expression for the specified gene
            df_gene = df[["celltype", "stage", gene]].copy()
            #groupby- groups cell types together and then groups by stage for
            #calculating below metrics
            #agg- new row with the row-wise mean and row-wise % cells expressing
            #pct_expr- where x > 0, True/1, else False/0; mean of this will be 
            # sum expressing/total cells, then *100
            #saving this as a separate df_summary
            df_summary = df_gene.groupby(["celltype", "stage"], observed=True).agg(
                mean_expr=(gene, "mean"), 
                pct_expr=(gene, lambda x: (x > 0).mean() * 100)
            ).reset_index()

            #now get the max % of cells expressing ligand for each cell type
            #across stages (for threshold)
            df_max_exp = df_summary.groupby(["celltype"], observed=True).agg(
                max_pct_exp=("pct_expr", "max")).reset_index()
            
            df_celltypes_keep = df_max_exp.copy(deep=True)[
                df_max_exp["max_pct_exp"] >= exp_threshold
            ]

            df_threshold = df_summary.copy(deep=True)[df_summary["celltype"].isin(list(df_celltypes_keep["celltype"]))]

            mean_expr_list.append(df_threshold["mean_expr"])
            pct_expr_list.append(df_threshold["pct_expr"])

    else:
        for gene in Hh_source_list:
            df_gene = df[["celltype", "stage", gene]].copy()
            df_summary = df_gene.groupby(["celltype", "stage"], observed=True).agg(
                mean_expr=(gene, "mean"),
                pct_expr=(gene, lambda x: (x > 0).mean() * 100)
            ).reset_index()
            mean_expr_list.append(df_summary["mean_expr"])
            pct_expr_list.append(df_summary["pct_expr"])

    mean_expr_all = pd.concat(mean_expr_list)
    pct_expr_all = pd.concat(pct_expr_list)
    mean_expr_min = mean_expr_all.min()
    mean_expr_max = mean_expr_all.max()
    pct_expr_min = pct_expr_all.min()
    pct_expr_max = pct_expr_all.max()

    #############
    # This code block actually creates the plots and subsets the data for plotting based 
    # on the specified threshold.
    #############
    # Set up subplots
    n_genes = len(Hh_source_list)
    if exp_threshold is not None:
        width_ = 7.5
    else:
        width_ = 4
    ##########################
    fig, axes = plt.subplots(ncols=n_genes, figsize=(width_*n_genes, fig_height), sharey=False)
    ##########################

    for ax, gene in zip(axes, Hh_source_list):
        dotplot_name = "dotplot_Hh_sources.svg"
        # this section computes the mean exp and % cells exp,
        # similarly to in the previous code block
        df_gene = df[["celltype", "stage", gene]].copy()
        df_summary = df_gene.groupby(["celltype", "stage"], observed=True).agg(
            mean_expr=(gene, "mean"),
            pct_expr=(gene, lambda x: (x > 0).mean() * 100)
        ).reset_index().copy(deep=True)

        if exp_threshold is not None:
            dotplot_name = f"dotplot_Hh_sources_has_exp_{exp_threshold}.svg"

            #now get the max % of cells expressing ligand for each cell type
            #across stages (for threshold)
            df_max_exp = df_summary.groupby(["celltype"], observed=True).agg(
                max_pct_exp=("pct_expr", "max")).reset_index()
            
            df_celltypes_keep = df_max_exp.copy(deep=True)[
                df_max_exp["max_pct_exp"] >= exp_threshold
            ]
            df_threshold = df_summary.copy(deep=True)[df_summary["celltype"].isin(list(df_celltypes_keep["celltype"]))]
            data_df = df_threshold.copy(deep=True)


        else:
            data_df = df_summary.copy(deep=True)

        hue_norm = mc.Normalize(vmin=mean_expr_min, vmax=mean_expr_max)
        size_norm = mc.Normalize(vmin=pct_expr_min, vmax=pct_expr_max)

        plotted_celltypes = sorted(data_df["celltype"].unique(), reverse=True)
        data_df["celltype"] = pd.Categorical(data_df["celltype"], categories=plotted_celltypes, ordered=True)

        # this block sets up the labeling of the fraction of cells expressing on the plot
        if annotate_frac_exp:
            # Compute number of expressing cells and total cells for each (celltype, stage)
            df_gene["is_expressing"] = df_gene[gene] > 0
            counts = df_gene.groupby(["celltype", "stage"], observed=True).agg(
            n_expressing=("is_expressing", "sum"),
            n_total=("is_expressing", "count")
            ).reset_index()
            # Merge counts with data_df to align order
            data_df = pd.merge(
            data_df,
            counts,
            on=["celltype", "stage"],
            how="left"
            )
            data_df["celltype"] = pd.Categorical(data_df["celltype"], categories=plotted_celltypes, ordered=True)

            # First: get the order of celltypes (same as what’s used in the plot)
            celltype_order = data_df["celltype"].cat.categories.tolist()

            # Map celltype string → y-axis index
            celltype_to_y = {ct: i for i, ct in enumerate(celltype_order)}

            # Annotate with vertical position based on numeric index
            for _, row in data_df.iterrows():
                y_pos = celltype_to_y[row["celltype"]] + 0.2
                ax.text(
                    row["stage"],
                    y_pos,
                    f"$\dfrac{{{int(row['n_expressing'])}}}{{{int(row['n_total'])}}}$",
                    ha="center",
                    va="bottom",
                    fontsize=6,
                    color="black"
                )
        size_min = 5
        size_max = 100
        size_range = size_max - size_min
        sns.scatterplot(
            data=data_df,
            x="stage",
            y="celltype",
            size="pct_expr",
            hue="mean_expr",
            sizes=(size_min, size_max),
            palette="Blues",
            edgecolor="black",
            linewidth=0.5,
            hue_norm=hue_norm,
            size_norm=size_norm,
            ax=ax,
            legend=False  # disable auto-legend
        )
        celltype_labels = data_df["celltype"].copy(deep=True)
        celltype_labels = celltype_labels.cat.rename_categories(
            [label.replace(" ", "\n") if " " in label else label for label in celltype_labels.cat.categories]
            )
        celltype_labels = celltype_labels.cat.rename_categories(
            [label.replace("/", "/\n") if "/" in label else label for label in celltype_labels.cat.categories]
            )
        
        xtick_vals = np.arange(len(df["stage"].cat.categories))
        ytick_vals = np.arange(len(celltype_labels.cat.categories))
        ax.set_title(f"{gene}", fontsize=8)
        ax.set_xlabel("Stage")
        ax.set_xticks(xtick_vals)
        ax.set_xlim(xtick_vals[0]-0.5, xtick_vals[-1]+0.5)
        ax.set_xticklabels(df["stage"].cat.categories),

        ax.set_yticks(ytick_vals)
        ax.set_yticklabels(celltype_labels.cat.categories)
        ax.set_ylabel("Cell Type")
        ##########################
        # add extra room at top #for plotting full dataset use 0.4, else 0.2
        ax.set_ylim(-0.2, len(plotted_celltypes) - 0.2)
        ##########################
    # Shared colorbar (for hue)
    sm = plt.cm.ScalarMappable(cmap="Blues", norm=hue_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation="vertical", shrink=0.6, pad=0.01)
    cbar.set_label("Mean Expression")
    pos = cbar.ax.get_position()
    cbar.ax.set_position([pos.x0, pos.y0 - 0.15, pos.width, pos.height])

    # Shared size legend (for % expressing)
    # Make dummy scatter handles with scaled sizes
    size_legend_vals = [pct_expr_min, (pct_expr_min + pct_expr_max) / 2, pct_expr_max]
    size_handles = [
        plt.scatter([], [], s=size_norm(val) * size_range + size_min,
                    color="gray", label=f"{int(val)}%")
        for val in size_legend_vals
    ]

    axes[-1].legend(
        handles=size_handles,
        fontsize=6,
        title="% Cells Expressing",
        title_fontsize=6,
        bbox_to_anchor=(1.05, 1.05),
        loc="upper left",
        borderaxespad=0.02,
        borderpad=0.75,
        labelspacing=1.1,         # ← add more vertical space between dots
        handletextpad=0.5,        # ← space between dot and label
        handlelength=1.3         # ← increase reserved space for dot
    )
    ##########################
    fig.subplots_adjust(wspace=0.35, hspace=0.35, right=0.77) #for full dataset, use wspace=0.4, else 0.35
    ##########################
    plt.savefig(fig_output_path+dotplot_name, bbox_inches='tight')
    # plt.show()
    plt.close()


def plot_stacked_violins(
        adata_dict: dict,
        Hh_source_list: list,
        fig_output_path: str
):
    """Creates a stacked violin plot of Hh_source
    ligand expression across stages for each cell 
    type. Makes one figure for each specified ligand.

    Parameters
    ----------
    adata_dict
        a dictionary of anndata objects, each containing
        the filtered Pijuan data for the given stage
        defined in the key

    Hh_source_list
        a list of strings defining the Hh ligands to plot

    fig_output_path
        a string defining the desired output directory for
        the plot
        ***needs to end in /***"""

    celltype_name = "celltype_PijuanSala2019"

    # calculate global medians across Hh sources
    # for violin colorbar scaling
    all_medians = []
    for gene in Hh_source_list:
        for stage, adata in adata_dict.items():
            labels = adata.obs[celltype_name].astype(str).values
            expr = adata.raw[:, gene].X.toarray().flatten()

            for ct in np.unique(labels):
                ct_mask = labels == ct
                all_medians.append(np.median(expr[ct_mask]))

    global_min = np.min(all_medians)
    global_max = np.max(all_medians)

    # get all unique celltypes across all time points
    # for plotting to have consistent y-axis labels/rows
    # at all time points
    for gene in Hh_source_list:
        all_celltypes = set()
        for stage, adata in adata_dict.items():
            labels = adata.obs[celltype_name].astype(str).values
            expr = adata.raw[:, gene].X.toarray().flatten()

            all_celltypes.update(labels)

        plotted_celltypes = sorted(all_celltypes)


        #set up the plot
        n_stages = len(list(adata_dict.keys()))
        fig, axes = plt.subplots(1, n_stages, figsize=(n_stages*0.75, 8)) #n_stages*0.75, *1.5 for wide
        if n_stages == 1:
            axes = [axes]

        for i, (ax, (stage, adata_stage)) in enumerate(zip(axes, adata_dict.items())):
            adata_stage = adata_stage.copy()

            # need to see if all cell types are included at this time point
            #if not, create dummmy data for that cell type so the y-axis
            #label and row still show up on the plot
            current_celltypes = set(adata_stage.obs[celltype_name].unique())
            missing_celltypes = set(plotted_celltypes) - current_celltypes

            # Add one dummy cell per missing cell type
            if missing_celltypes:
                n_genes = adata_stage.raw.shape[1] if adata_stage.raw is not None else adata_stage.shape[1]
                #use nans so the violin doesn't actually show up on plot
                dummy_X = np.full((len(missing_celltypes), n_genes), np.nan)
                dummy_obs = pd.DataFrame({
                    celltype_name: list(missing_celltypes)
                }, index=[f"dummy_{ct}" for ct in missing_celltypes])

                # Create dummy AnnData
                dummy_adata = sc.AnnData(
                    X=dummy_X,
                    obs=dummy_obs,
                    var=adata_stage.raw.var if adata_stage.raw is not None else adata_stage.var
                )

                # Combine with original
                adata_stage = ad.concat([adata_stage, dummy_adata])

                # Restore raw if needed
                if adata_stage.raw is not None:
                    adata_stage.raw = adata_stage[:, adata_stage.raw.var_names]

            adata_stage.obs[celltype_name] = pd.Categorical(
                adata_stage.obs[celltype_name].copy(),
                categories=plotted_celltypes,
                ordered=True
            )
            
            vp = sc.pl.stacked_violin(
                adata_stage,
                var_names=[gene],
                groupby=celltype_name,
                stripplot=False,
                ax=ax,
                show=False,
                return_fig=True,
                vmin=global_min,
                vmax=global_max,
            )
            axes_dict = vp.get_axes()
            axes_dict["mainplot_ax"].set_xticklabels([stage], rotation=0)
            axes_dict["mainplot_ax"].tick_params(axis='y', labelsize=6)
            axes_dict["mainplot_ax"].tick_params(axis='x', labelsize=6)
            axes_dict["mainplot_ax"].set_ylabel("Cell type")

            if i > 0:
                axes_dict["mainplot_ax"].yaxis.set_visible(False)
                axes_dict["mainplot_ax"].spines['left'].set_visible(False)

            # Remove extra colorbar
            bar_ax = fig.axes[-1]
            if bar_ax not in axes:
                bar_ax.remove()

        # Step 4: Shared colorbar
        sm = plt.cm.ScalarMappable(cmap="Blues", norm=mc.Normalize(vmin=global_min, vmax=global_max))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, orientation="vertical", shrink=0.3, pad=0.01)
        cbar.set_label("Median Expression")
        pos = cbar.ax.get_position()
        cbar.ax.set_position([pos.x0 - 0.1, pos.y0, pos.width, pos.height]) #x + 0.075, y - 0.125
        
        plot_name = f"violins_by_celltype_{gene}.svg"
        fig.suptitle(gene, fontsize=7, y=0.85, x=0.37)
        fig.supxlabel("Stage", fontsize=6, y=0.07, x=0.37)
        fig.subplots_adjust(wspace=-0.85, hspace=0.0)
        plt.savefig(fig_output_path+plot_name, bbox_inches='tight')
        plt.close()
