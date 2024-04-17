import marimo

__generated_with = "0.4.0"
app = marimo.App()


@app.cell
def __(
    currcat,
    genes_in_category,
    level,
    mo,
    neur_fig,
    neuron,
    path_fig,
    pathway,
):
    tab1 = mo.vstack(
        [
            mo.hstack([level, neuron, pathway], justify="start"),
            mo.ui.plotly(neur_fig),
            mo.ui.table(genes_in_category[pathway.value]) if pathway.value else "",
        ]
    )

    tab2 = mo.vstack([level,currcat, mo.ui.plotly(path_fig)])

    tabs = mo.ui.tabs(
        {
            "Neuron Analysis": tab1,
            "Pathway across neurons": tab2,
        }
    )
    tabs
    return tab1, tab2, tabs


@app.cell
def __():
    import marimo as mo
    import pandas as pd
    import re
    import plotly.express as px
    return mo, pd, px, re


@app.cell
def __(pd):
    ## DATA IMPORT
    genes = pd.read_csv("data/genes.csv").set_index("WormBase ID")
    #print(genes)
    cats = pd.read_csv("data/four_level_categories.csv")
    L1_genes = pd.read_csv("data/level_one_genes.csv")
    L2_genes = pd.read_csv("data/level_two_genes.csv")
    L3_genes = pd.read_csv("data/level_three_genes.csv")
    L4_genes = pd.read_csv("data/level_four_genes.csv")

    # lookup function:
    # genes.loc[genes['Gene Name'] == 'pfk-1.1', 'WormBase ID'].iloc[0]

    category_genes = {
        "LEVEL 1": L1_genes,
        "LEVEL 2": L2_genes,
        "LEVEL 3": L3_genes,
        "LEVEL 4": L4_genes,
    }

    reads = pd.read_csv("data/NormToMed_Bulk-CenGen_111521.csv")
    for _col in reads.columns[1:]:
        reads[_col] = pd.to_numeric(reads[_col], errors="coerce")

    # bwm = pd.read_csv("data/NonBulk_BWM.csv")
    # hyp = pd.read_csv("data/NonBulk_HYP.csv")

    # reads.set_index("Gene")

    # _filtered = reads.index.intersection(_genes)]
    return (
        L1_genes,
        L2_genes,
        L3_genes,
        L4_genes,
        category_genes,
        cats,
        genes,
        reads,
    )


@app.cell
def __(pd, re, reads):
    ## Import the normalized data and average for each neuron type.
    _curr_neur = None
    _tmp = pd.DataFrame()
    avg = pd.DataFrame()


    ## Loops through the DataFrame and checks the column names, splitting at "r".
    #  It averages the columns in a dataframe, then adds the result into the "avg" dataframe
    for _col in reads.columns[1:]:
        _neur = re.split("[r]", _col)[0]
        if _neur != _curr_neur and _curr_neur != None:
            avg.insert(len(avg.columns), _curr_neur, _tmp.mean(axis=1))
            _tmp = pd.DataFrame()
            _curr_neur = _neur
            _tmp.insert(len(_tmp.columns), len(_tmp.columns), reads[_col])

        elif _neur == _curr_neur:
            _tmp.insert(len(_tmp.columns), len(_tmp.columns), reads[_col])
        else:
            _curr_neur = _neur
            _tmp.insert(len(_tmp.columns), len(_tmp.columns), reads[_col])

    avg.insert(0, "Gene", reads["Gene"])
    avg.set_index("Gene", inplace=True)
    return avg,


@app.cell
def __(avg, category_genes, genes, level, mo, neuron, pd, px):
    ## Iterate through the categories at the chosen level
    #  for each category, look up the RNA value for each gene
    #  then place this in a new dataframe
    _effects = {}
    genes_in_category = {}

    for _col in category_genes[level.value].columns:
        _genes = category_genes[level.value][_col].dropna()
        _filtered = avg[neuron.value].loc[
            avg[neuron.value].index.intersection(_genes)
        ]
        _effects[_col] = _filtered.dropna().mean()
        _names = genes.loc[genes.index.intersection(_filtered.index)]
        _combined = pd.concat([_names, _filtered], axis=1)
        genes_in_category[_col] = _combined

    _effects_sorted = dict(sorted(_effects.items(), key=lambda item: item[1]))

    _ht = 400 if level.value == "LEVEL 1" else 1500

    neur_fig = px.scatter(
        x=_effects_sorted.values(), y=_effects_sorted.keys(), height=_ht
    )
    neur_fig.update_xaxes(type="log")
    neur_fig.update_yaxes(type="category")
    neur_fig.update_layout(yaxis={"dtick": 1})
    mo.output.clear()
    return genes_in_category, neur_fig


@app.cell
def __(avg, mo):
    level = mo.ui.dropdown(
        ["LEVEL 1", "LEVEL 2", "LEVEL 3", "LEVEL 4"], value="LEVEL 1"
    )
    neuron = mo.ui.dropdown(avg.columns[1:], value="ASER")
    return level, neuron


@app.cell
def __(category_genes, level, mo):
    pathway = mo.ui.dropdown(category_genes[level.value].columns)
    return pathway,


@app.cell
def __(category_genes, level, mo):
    currcat = mo.ui.dropdown(
        category_genes[level.value].columns,
        value=category_genes[level.value].columns[0],
    )
    return currcat,


@app.cell
def __(avg, category_genes, currcat, level, mo, px):
    _effects = {}

    _genes = category_genes[level.value][currcat.value].dropna()
    _filtered = avg.loc[avg.index.intersection(_genes)]

    for _col in _filtered.columns:
        _effects[_col] = _filtered[_col].mean()

    _effects_sorted = dict(sorted(_effects.items(), key=lambda item: item[1]))

    path_fig = px.scatter(
        x=_effects_sorted.values(), y=_effects_sorted.keys(), height=700
    )
    path_fig.update_xaxes(type="log")
    path_fig.update_yaxes(type="category")
    path_fig.update_layout(yaxis={"dtick": 1})
    mo.output.clear()
    return path_fig,


if __name__ == "__main__":
    app.run()
