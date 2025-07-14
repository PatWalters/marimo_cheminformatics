import marimo

__generated_with = "0.14.10"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo, txta):
    txt = mo.ui.text("")
    txta
    return (txt,)


@app.cell
def _(txt):
    txt.value
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
