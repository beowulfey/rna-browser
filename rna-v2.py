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