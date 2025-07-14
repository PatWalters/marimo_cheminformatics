import marimo

__generated_with = "0.14.10"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## REOS in a Marimo Notebook
    This notebook provides a quick overview of how the [Rapid Elimination of Swill (REOS)](https://practicalcheminformatics.blogspot.com/2018/08/filtering-chemical-libraries.html) filters can be run in a Marimo notebook.  We'll take a look at how the [useful_rdkit_utils](https://patwalters.github.io/Useful-RDKit-Utils/) library and Marimo's radio button capability can be used to create a quick interactive viewer.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""First, we'll import the necessary Python libraries.""")
    return


@app.cell
def _():
    from rdkit import Chem
    import useful_rdkit_utils as uru
    import pandas as pd
    import marimo as mo
    return Chem, mo, pd, uru


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Define a quick function to return the SMILES for the largest fragment in a molecule.""")
    return


@app.cell
def _(Chem, uru):
    def strip_salts(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol =  uru.get_largest_fragment(mol)
            return Chem.MolToSmiles(mol)
        else:
            return None
    return (strip_salts,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ##### Read and process a SMILES file.      
    1. Read the SMILES file into a Pandas dataframe.     
    2. Use the function we defined about to remove counterions, waters of hydration, etc.     
    3. Use the [smi_to_base64_image](https://useful-rdkit utils.readthedocs.io/en/latest/misc_utils.html#useful_rdkit_utils.misc_utils.smi_to_base64_image) method in `useful_rdkit_utils` to add a structure image to the dataframe.  Marimo will show an expanded version of the image as tooltip when the mouse is hovered over the structure image.
    """
    )
    return


@app.cell
def _(pd, strip_salts, uru):
    df = pd.read_csv("test.smi",sep=" ",names=["SMILES","Name"])
    df.SMILES = df.SMILES.apply(strip_salts)
    df.insert(0,'image',df.SMILES.apply(uru.smi_to_base64_image,target="altair"))
    return (df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Instantiate a REOS object.""")
    return


@app.cell
def _(uru):
    reos = uru.REOS()
    return (reos,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Run the filters on `df`.  Note that this applies the REOS SMARTS patterns to each of the SMILES in the dataframe.""")
    return


@app.cell
def _(df, reos):
    reos_df = reos.pandas_smiles(df.SMILES)
    df['reos'] = reos_df.description
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""The useful_rdkit_utils library has a method called [value_counts_df](https://useful-rdkit-utils.readthedocs.io/en/latest/pandas.html#useful_rdkit_utils.pandas_utils.value_counts_df) that creates a new dataframe with the counts of unique entries in another dataframe.  Under the hood, this is just running the Pandas `value_counts` method and putting the results into a dataframe.""")
    return


@app.cell
def _(df, uru):
    reos_summary_df = uru.value_counts_df(df,"reos")
    reos_summary_df
    return (reos_summary_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""With `reos_summary_df` in hand, we can use Marimo's `radio` component to create an interactive viewer.  First, we'll pass the list of REOS filtlers that were triggered by the dataset and use this to create a radio button panel called `filter`.  We can get the selected radio button with `filter.value`.""")
    return


@app.cell
def _(mo, reos_summary_df):
    filter = mo.ui.radio(reos_summary_df.reos)
    return (filter,)


@app.cell
def _(df, filter, mo):
    mo.hstack([filter,df.query(f"reos == '{filter.value}'")])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
