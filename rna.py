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
def __(cell, cell_fig, mo, path_fig, pathway, pathway_genes_converted):
    ## Some UI stuff

    cell_tab = mo.vstack(
        [
            cell,
            cell_fig
        ]
    )

    path_tab = mo.vstack(
        [
            pathway,
            path_fig
        ]
    )

    info_tab = mo.vstack(
        [
            pathway,
            mo.ui.table(
                pathway_genes_converted,
                pagination=True,
                page_size=10,
            ),
        ]
    )

    tabs = mo.ui.tabs(
        {
            "Pathway Enrichment per Cell": cell_tab,
            "Cell Enrichment per Pathway": path_tab,
            "Pathway Gene Lists": info_tab,
        }
    )

    tabs
    return cell_tab, info_tab, path_tab, tabs


@app.cell
def __(mo):
    ## UI elements
    source = mo.ui.dropdown(["Classic CeNGEN", "Bulk CeNGEN"], value="Bulk CeNGEN")
    level = mo.ui.dropdown(["LEVEL 1", "LEVEL 2", "LEVEL 3", "LEVEL 4"], value="LEVEL 1")
    # neuron = mo.ui.dropdown(reads.columns[1:], value="AIY")
    return level, source


@app.cell
def __(category_genes, level, mo, reads):
    ## UI Elements that depend on other elements

    pathway = mo.ui.dropdown(category_genes[level.value].columns, value=category_genes[level.value].columns[0])
    cell = mo.ui.dropdown(reads.columns.to_list(), value=reads.columns.to_list()[0])

    currcat = mo.ui.dropdown(
        category_genes[level.value].columns,
        value=category_genes[level.value].columns[0],
    )
    return cell, currcat, pathway


@app.cell
def __(category_genes, genes, level, pathway):
    ## UI Elements: third tier dependence

    pathway_genes_converted = genes.loc[genes.index.intersection(category_genes[level.value][pathway.value].dropna())]
    #pathway_genes_converted.rename("Gene", inplace=True)
    #pathway_genes_converted
    #mo.output.clear()
    return pathway_genes_converted,


@app.cell
def __(level, pd, re, source):
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
    url_genes = "https://raw.githubusercontent.com/beowulfey/rna-browser/main/data/wormbase_genes.csv"


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

    if source.value == "Bulk CeNGEN":
        _curr_neur = None
        _tmp = pd.DataFrame()
        _avg = pd.DataFrame()
        for _col in reads.columns:
            _neur = re.split("[r]", _col)[0]
            if _neur != _curr_neur and _curr_neur != None:
                _avg.insert(len(_avg.columns), _curr_neur, _tmp.mean(axis=1))
                _tmp = pd.DataFrame()
                _curr_neur = _neur
                _tmp.insert(len(_tmp.columns), len(_tmp.columns), reads[_col])
        
            elif _neur == _curr_neur:
                _tmp.insert(len(_tmp.columns), len(_tmp.columns), reads[_col])
            else:
                _curr_neur = _neur
                _tmp.insert(len(_tmp.columns), len(_tmp.columns), reads[_col])
        reads = _avg


    for _col in reads.columns:
        reads[_col] = pd.to_numeric(reads[_col], errors="coerce")


    genes = pd.read_csv(url_genes)
    genes.set_index("Wormbase_ID", inplace=True)

    cells = reads.iloc[0].index.tolist()

    pathway_genes = {}
    for col in category_genes[level.value].columns:
        pathway_genes[col] = category_genes[level.value][col].dropna().tolist()
    return (
        L1_genes,
        L2_genes,
        L3_genes,
        L4_genes,
        category_genes,
        cats,
        cells,
        col,
        genes,
        pathway_genes,
        reads,
        url_1,
        url_2,
        url_3,
        url_4,
        url_cats,
        url_genes,
        url_sc,
    )


@app.cell
def __(cell, np, pathway_genes, pathway_genes_converted, reads):
    per_cell = {}
    cell_sum = {}
    path_sum = {}
    meds = {}

    for _gene in pathway_genes_converted.index:
        if _gene in reads.index:
            _med = np.median([i for i in reads.loc[_gene].to_list() if i >= 0])
            meds[_gene] = _med
            for _cell in reads.loc[_gene].index:
                if not _cell in per_cell.keys():
                    per_cell[_cell] = {}
                if _med > 0:
                    per_cell[_cell][_gene] = (reads.loc[_gene].loc[_cell]) / _med
                else:
                    per_cell[_cell][_gene] = 0
    for _cell in per_cell:
        cell_sum[_cell] = np.median([i for i in list(per_cell[_cell].values()) if i >= 0])

    for (_path,_genes) in pathway_genes.items():
        _tmp_path = []
        for _gene in _genes:
            if _gene in reads.index:
                _tmp_path.append(reads[cell.value].loc[_gene])
        path_sum[_path] = np.median(_tmp_path)

    return cell_sum, meds, path_sum, per_cell


@app.cell
def __(cell, cell_sum, mo, path_sum, pathway, px):
    _sorted_cells = {k: v for k, v in sorted(cell_sum.items(), key=lambda item: item[1])}
    _filtered_cells = {k: v for k, v in _sorted_cells.items() if v > 0}
    path_fig = px.scatter(x=_filtered_cells.values(), y=_filtered_cells.keys(), height=len(_filtered_cells.values()) * 20, title=pathway.value)
    path_fig.update_xaxes(type="log")
    path_fig.update_yaxes(type="category")
    path_fig.update_layout(yaxis={"dtick": 1})

    _sorted_paths = {k: v for k, v in sorted(path_sum.items(), key=lambda item: item[1])}
    _filtered_paths = {k: v for k, v in _sorted_paths.items() if v > 0}
    cell_fig = px.scatter(x=_filtered_paths.values(), y=_filtered_paths.keys(), height=len(_filtered_paths.values()) * 20, title=cell.value)
    cell_fig.update_xaxes(type="log")
    cell_fig.update_yaxes(type="category")
    cell_fig.update_layout(yaxis={"dtick": 1})

    mo.output.clear()
    return cell_fig, path_fig


@app.cell
def __():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import plotly.express as px
    import re

    # making sure plots & clusters are reproducible
    np.random.seed(42)
    return mo, np, pd, px, re


if __name__ == "__main__":
    app.run()
