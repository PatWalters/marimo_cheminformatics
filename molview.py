import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import useful_rdkit_utils as uru
    from rdkit import Chem
    return Chem, uru


@app.cell
def _(uru):
    uru.rd_enable_svg()
    return


@app.cell
def _(Chem):
    mol_list = [Chem.MolFromSmiles("CCCC")]
    return (mol_list,)


@app.cell
def _(mol_list):
    mol_list[0]
    return


@app.cell
def _(mol_list, uru):
    uru.mol_to_3D_view(mol_list)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
