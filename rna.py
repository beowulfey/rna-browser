import marimo

__generated_with = "0.5.2"
app = marimo.App(width="medium")


@app.cell
def __(level, mo, source):
    mo.vstack(
        [
            mo.md(
                """
        ## Neuronal metabolism browser
        Uses CenGEN data to compare metabolic profiles of different neurons.   
        Data adapted from https://cengen.shinyapps.io/CengenApp/ and has data for all neurons (threshold 2).   
        """
            ),
            mo.hstack(["Data source: ", source], justify="start"),
            mo.hstack(["Pathway resolution: ", level], justify="start"),
        ]
    )
    return


@app.cell
def __(mo, pathway, pathway_genes_converted):
    ## Some UI stuff

    data_tab = mo.vstack(
        [
            pathway,
        ]
    )

    info_tab = mo.vstack(
        [
            pathway,
            mo.ui.table(
                pathway_genes_converted.to_frame(),
                pagination=True,
                page_size=10,
            ),
        ]
    )

    tabs = mo.ui.tabs(
        {
            "Data Explorer": data_tab,
            "Pathway Gene Lists": info_tab,
        }
    )

    tabs
    return data_tab, info_tab, tabs


@app.cell
def __(mo):
    ## UI elements
    source = mo.ui.dropdown(["Classic CeNGEN", "Bulk CeNGEN"], value="Bulk CeNGEN")
    level = mo.ui.dropdown(["LEVEL 1", "LEVEL 2", "LEVEL 3", "LEVEL 4"], value="LEVEL 1")
    # neuron = mo.ui.dropdown(reads.columns[1:], value="AIY")
    return level, source


@app.cell
def __(category_genes, level, mo):
    ## UI Elements that depend on other elements

    pathway = mo.ui.dropdown(category_genes[level.value].columns, value=category_genes[level.value].columns[0])

    currcat = mo.ui.dropdown(
        category_genes[level.value].columns,
        value=category_genes[level.value].columns[0],
    )
    return currcat, pathway


@app.cell
def __(category_genes, genes, level, mo, pathway):
    ## UI Elements: third tier dependence

    pathway_genes_converted = genes.loc[genes.index.intersection(category_genes[level.value][pathway.value].dropna())]
    pathway_genes_converted.rename("Gene", inplace=True)

    mo.output.clear()
    return pathway_genes_converted,


@app.cell
def __(pd, source):
    ## DATA IMPORT
    #  Imports the categories (lists of genes in each category) at four levels of detail. Imports the CenGen Reads and pulls the gene names
    #  out into a separate dateframe.

    ## CENGEN DATASET
    if source.value == "Bulk CeNGEN":
        url_sc = "https://raw.githubusercontent.com/beowulfey/rna-browser/main/data/Bulk-CenGen_111521.csv"
    else:
        url_sc = "https://raw.githubusercontent.com/beowulfey/rna-browser/main/data/Liska_Single_Cell_TPM_threshold2.csv"

    ## WALHOUT LAB ANNOTATIONS
    url_cats = "https://raw.githubusercontent.com/beowulfey/rna-browser/main/data/four_level_categories.csv"
    url_1 = "https://raw.githubusercontent.com/beowulfey/rna-browser/main/data/level_one_genes.csv"
    url_2 = "https://raw.githubusercontent.com/beowulfey/rna-browser/main/data/level_two_genes.csv"
    url_3 = "https://raw.githubusercontent.com/beowulfey/rna-browser/main/data/level_three_genes.csv"
    url_4 = "https://raw.githubusercontent.com/beowulfey/rna-browser/main/data/level_four_genes.csv"

    cats = pd.read_csv(url_cats)
    L1_genes = pd.read_csv(url_1)
    L2_genes = pd.read_csv(url_2)
    L3_genes = pd.read_csv(url_3)
    L4_genes = pd.read_csv(url_4)

    category_genes = {
        "LEVEL 1": L1_genes,
        "LEVEL 2": L2_genes,
        "LEVEL 3": L3_genes,
        "LEVEL 4": L4_genes,
    }

    reads = pd.read_csv(url_sc)
    reads.set_index("Wormbase_ID", inplace=True)
    for _col in reads.columns[:]:
        reads[_col] = pd.to_numeric(reads[_col], errors="coerce")

    genes = reads["Gene_Name"]
    reads.drop(["Sequence_Name", "Gene_Name"], axis=1, inplace=True)

    cells = reads.iloc[0].index.tolist()
    return (
        L1_genes,
        L2_genes,
        L3_genes,
        L4_genes,
        category_genes,
        cats,
        cells,
        genes,
        reads,
        url_1,
        url_2,
        url_3,
        url_4,
        url_cats,
        url_sc,
    )


@app.cell
def __(np, pathway_genes_converted, reads):
    per_cell = {}
    cell_sum = {}

    for gene in pathway_genes_converted.index:
        med = np.median([i for i in reads.loc[gene].to_list() if i > 0])

        for cell in reads.loc[gene].index:
            # if cell == 'ASER':
            # print(gene, pathway_genes_converted.loc[gene], cell, reads.loc[gene].loc[cell],reads.loc[gene].loc[cell]/med)
            if not cell in per_cell.keys():
                per_cell[cell] = {}
            if med > 0:
                per_cell[cell][gene] = (reads.loc[gene].loc[cell]) / med
            else:
                per_cell[cell][gene] = 0

    for cell in per_cell:
        cell_sum[cell] = np.median([i for i in list(per_cell[cell].values()) if i > 0])
        # cell_sum[cell] = sum(list(per_cell[cell].values()))/len(list(per_cell[cell].values()))
    return cell, cell_sum, gene, med, per_cell


@app.cell
def __(cell_sum, px):
    _sorted = {k: v for k, v in sorted(cell_sum.items(), key=lambda item: item[1])}
    _filtered = {k: v for k, v in _sorted.items() if v > 0}
    neur_fig = px.scatter(x=_filtered.values(), y=_filtered.keys(), height=len(_filtered.values()) * 20)
    neur_fig.update_xaxes(type="log")
    neur_fig.update_yaxes(type="category")
    neur_fig.update_layout(yaxis={"dtick": 1})
    return neur_fig,


@app.cell
def __():
    return


@app.cell
def __():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import plotly.express as px

    # making sure plots & clusters are reproducible
    np.random.seed(42)
    return mo, np, pd, px


if __name__ == "__main__":
    app.run()
