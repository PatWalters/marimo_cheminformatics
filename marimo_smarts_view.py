import marimo

__generated_with = "0.14.10"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### A simple SMILES/SMARTS viewer with [marimo](https://marimo.io/)

    This notebook provides a quick introduction to some of the reactive capbilities in **marimo**.  We will create simple viewer for SMILES and SMARTS.  When the user types in the text boxes below, the chemical structure associated with the SMILES and the highlighting based on the SMARTS witll automatically update.  

    We begin by importing the necessary python libraries.
    """
    )
    return


@app.cell
def _():
    from rdkit import Chem
    import useful_rdkit_utils as uru
    import marimo as mo
    return Chem, mo, uru


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""By default, **marimo** won't display RDKit Molecule objects as chemical structures.  However, this function from [useful_rdkit_utils](https://github.com/PatWalters/useful_rdkit_utils) will enable structure display.""")
    return


@app.cell
def _(uru):
    uru.rd_enable_svg()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Define a function to display the chemical structure represented by a SMILES and highlight the atoms represented by the SMARTS pattern in the second function argument.""")
    return


@app.cell
def _(Chem):
    def view_mol(smiles_val, smarts_val):
        mol = Chem.MolFromSmiles(smiles_val)
        if mol is not None:
            pat = Chem.MolFromSmarts(smarts_val)
            if pat:
                match = mol.GetSubstructMatch(pat)
                if len(match) == 0:
                    print("No match")
        return mol
    return (view_mol,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Define the text fields we use for SMILES and SMARTS input.  Here we simply define two [marimo.ui.text](https://docs.marimo.io/api/inputs/text/) fields. We can get the values entered into the text fields using `smiles.value` and `smarts.value`.  Note that `debounce` is set to `True`, so the associated value will update immediately as the user types.  If we set `debounce` to `False`, the values won't be updated until the user presses Enter or the text field loses focus.""")
    return


@app.cell
def _(mo):
    smiles = mo.ui.text(label="SMILES",debounce=False)
    smarts = mo.ui.text(label="SMARTS",debounce=False)
    return smarts, smiles


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    In the cell below we define the viewer.  We use the marimo [vstack]() function to create a vertical "stack" of three components.    
    - The SMILES text field   
    - The SMARTS text field   
    - The molecule returned by the `view_mol` function above

    Give it a try.  Note that the structure display updates as you type.       
    If you're look for ideas, try copying and pasting these: 
    ```
    SMILES: CN1CCN(CC1)Cc1ccc(cc1)C(=O)Nc1ccc(c(c1)Nc1nccc(n1)c1cccnc1)C
    SMARTS: c1ccccn1
    ```
    """
    )
    return


@app.cell
def _(mo, smarts, smiles, view_mol):
    mo.vstack([smiles,smarts,view_mol(smiles.value, smarts.value)])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
