import marimo

__generated_with = "0.4.0"
app = marimo.App(width="medium")


@app.cell
def __(level, mo):
    mo.vstack(
        [
            mo.md(
                """
        ## Neuronal metabolism browser
        Uses CenGEN data to compare metabolic profiles of different neurons.   
        Data adapted from https://cengen.shinyapps.io/CengenApp/ and has data for all neurons (threshold 2).   
        """
            ),
            mo.hstack(["Pathway resolution: ", level], justify="start"),
        ]
    )
    return


@app.cell
def __(mo, pathway, pathway_genes_converted):
    ## Some UI stuff

    data_tab = "PLACEHOLDER"

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
def __():
    ## scRNA ANALYSIS 1
    #  This uses Meld to analyze the different pathways.

    # Lets try going straight from the raw data
    sample_info = [('GSM4040097', 'Pan-neural', '1'), ('GSM4040098', 'Pan-neural', '2'), ('GSM4040099', 'GABAergic neurons', '1'),
               ('GSM404100', 'nmr-1 expressing neurons', '1'), ('GSM404101', 'Glutamatergic neurons', '1'), ('GSM404102', 'GABAergic neurons', '2')]


    return sample_info,


@app.cell
def __():
    return


@app.cell
def __(mo, reads):
    ## UI elements

    level = mo.ui.dropdown(["LEVEL 1", "LEVEL 2", "LEVEL 3", "LEVEL 4"], value="LEVEL 1")
    neuron = mo.ui.dropdown(reads.columns[1:], value="AIY")
    return level, neuron


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
    pathway_genes_converted.rename(pathway.value, inplace=True)

    mo.output.clear()
    return pathway_genes_converted,


@app.cell(disabled=True)
def __(pd):
    ## DATA IMPORT V1
    #  Imports the categories (lists of genes in each category) at four levels of detail. Imports the CenGen Reads and pulls the gene names
    #  out into a separate dateframe.


    ## CENGEN DATASET
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
    for _col in reads.columns[4:]:
        reads[_col] = pd.to_numeric(reads[_col], errors="coerce")

    genes = reads["Gene_Name"]
    reads.drop(["Sequence_Name", "Gene_Name"], axis=1, inplace=True)
    return (
        L1_genes,
        L2_genes,
        L3_genes,
        L4_genes,
        category_genes,
        cats,
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
def __():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import scprep
    import meld
    import plotly.express as px

    # making sure plots & clusters are reproducible
    np.random.seed(42)


    return meld, mo, np, pd, px, scprep


if __name__ == "__main__":
    app.run()
