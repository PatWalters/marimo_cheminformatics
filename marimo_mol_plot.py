import marimo

__generated_with = "0.14.10"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    import useful_rdkit_utils as uru
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    import altair as alt
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    import base64
    import marimo as mo
    from rdkit.Chem.Draw import MolsToGridImage
    return (
        Chem,
        MolsToGridImage,
        PCA,
        TSNE,
        alt,
        base64,
        mo,
        pd,
        rdMolDraw2D,
        uru,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""To begin, we select a SMILES file. We can do that by either clicking on the button below or dragging a file onto the button.""")
    return


@app.cell
def _():
    filename = "test.smi"
    return (filename,)


@app.cell
def _(filename, pd):
    df = pd.read_csv(filename,sep=" ",names=["SMILES","Name"])
    return (df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Use the Rule of 5 calculator from useful_rdkit_utils to calculate properties. Once properties are calculated, add them to `df`.""")
    return


@app.cell
def _(df, uru):
    ro5 = uru.Ro5Calculator()
    df_ro5 = ro5.pandas_smiles(df.SMILES)
    for col in df_ro5.columns:
        df[col] = df_ro5[col]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Use the `Smi2Fp` capability in `useful_rdkit_utils` to calculate fingerprints.""")
    return


@app.cell
def _(df, uru):
    smi2fp = uru.Smi2Fp()
    df["fp"] = df["SMILES"].apply(smi2fp.get_np_counts)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Use tSNE to create a 2D projection of the fingerprint space.""")
    return


@app.cell
def _(PCA, TSNE, df):
    pca = PCA(n_components=50)
    pcs = pca.fit_transform(list(df["fp"]))
    tsne = TSNE(n_components=2, perplexity=30, max_iter=300)
    tsne_result = tsne.fit_transform(pcs)
    return (tsne_result,)


@app.cell
def _(df, tsne_result):
    df[['tsne_x','tsne_y']] = tsne_result.tolist()
    return


@app.cell
def _(Chem, base64, rdMolDraw2D):
    def mol_to_base64_image(mol: Chem.Mol) -> str:
        """
        Convert an RDKit molecule to a base64 encoded image string.

        Parameters:
        mol (Chem.Mol): The RDKit molecule to convert.

        Returns:
        str: The base64 encoded image string.
        """
        drawer = rdMolDraw2D.MolDraw2DCairo(300, 150)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        text = drawer.GetDrawingText()
        im_text64 = base64.b64encode(text).decode('utf8')
        return f"data:image/png;base64,{im_text64}"
    return (mol_to_base64_image,)


@app.cell
def _(Chem, mol_to_base64_image):
    def smi2image(smi):
        mol = Chem.MolFromSmiles(smi)
        return mol_to_base64_image(mol)
    return


@app.cell
def _(df, uru):
    df['image'] = df.SMILES.apply(uru.smi_to_base64_image,target="altair")
    return


@app.cell
def _(alt, df):
    alt.Chart(df).mark_circle(size=60).encode(
        x = alt.X('tsne_x',title='t-SNE X',
            axis=alt.Axis(
                titleFontSize=16,  # Font size for the x-axis title
                labelFontSize=12   # Font size for the x-axis tick labels (numbers/text)
            )),
        y = alt.Y('tsne_y',title='t-SNE Y',
            axis=alt.Axis(
                titleFontSize=16,  # Font size for the x-axis title
                labelFontSize=12   # Font size for the x-axis tick labels (numbers/text)
            )),
        color=alt.condition(
            alt.datum.LogP < 5,  # Condition: if the 'value' field is less than 5
            alt.value('blue'),    # If true, set color to blue
            alt.value('red')      # If false (i.e., value >= 5), set color to red
        ),
        tooltip=['image','Name']).interactive()
    return


@app.cell
def _(df, mo):
    skip_columns = ['SMILES','fp','tsne_x','tsne_y','image']
    show_cols = ['image'] + [x for x in df.columns if x not in skip_columns]
    mo_df = mo.ui.table(df[show_cols].round(2))
    mo_df
    return mo_df, show_cols


@app.cell
def _(mo_df):
    len(mo_df.value)
    return


@app.cell
def _(alt, df, mo):
    chart = mo.ui.altair_chart(alt.Chart(df).mark_circle(size=60).encode(
        x = alt.X('tsne_x',title='t-SNE X',
            axis=alt.Axis(
                titleFontSize=16,  # Font size for the x-axis title
                labelFontSize=12   # Font size for the x-axis tick labels (numbers/text)
            )),
        y = alt.Y('tsne_y',title='t-SNE Y',
            axis=alt.Axis(
                titleFontSize=16,  # Font size for the x-axis title
                labelFontSize=12   # Font size for the x-axis tick labels (numbers/text)
            )),
        color=alt.condition(
            alt.datum.LogP < 5,  # Condition: if the 'value' field is less than 5
            alt.value('blue'),    # If true, set color to blue
            alt.value('red')      # If false (i.e., value >= 5), set color to red
        ),
        tooltip=['image','Name']))
    return (chart,)


@app.cell
def _(chart, mo, show_cols):
    mo.vstack([chart,chart.value[show_cols].round(1)])
    return


@app.cell(hide_code=True)
def _(chart):
    chart.value.SMILES.values
    return


@app.cell
def _(Chem, MolsToGridImage):
    def display_smiles(smi_list):
        mol_list = [Chem.MolFromSmiles(x) for x in smi_list]
        if len(mol_list):
            return MolsToGridImage(mol_list, molsPerRow=8, subImgSize=(150, 150))
        else:
            return None
    return (display_smiles,)


@app.cell
def _(chart, display_smiles, mo):
    mo.vstack([chart,display_smiles(chart.value.SMILES.values[:8])])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
