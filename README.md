# marimo_cheminformatics
Practical Cheminformatics with Marimo

The notebooks in this repo show some examples of how to use marimo notebooks as a tool for cheminformatic analysis. Marimo is often referred to as a "better Jupyter notebook". I've found two primary advantages of marimo over Jupyter.
- Marimo is reactive, changes to a variable in one cell will be updated all other cells containing that variable.
- Marimo contains several native GUI components that can be integrated to easily generate interactive analyses. 

For more information on marimo, I'd recommend these tutorials and references.
- [Marimo's getting started docs](https://docs.marimo.io/getting_started/)
- [Marimo FAQ](https://docs.marimo.io/faq/)
- [Marimo article in Real Python](https://realpython.com/marimo-notebook/)


### Installing marimo
The marimo docs provide example of how to install marimo with [pip](https://docs.marimo.io/#__tabbed_1_1),
[uv](https://docs.marimo.io/#__tabbed_1_2), or [conda](https://docs.marimo.io/#__tabbed_1_3).

### Running the notebooks
The notebooks can be run with `marimo edit notebook.py`, where `notebook.py` is the 

### The current collection
This repo currently contains a few notebooks demonstrating how to use marimo for cheminformatics. I'll continue add notebooks. If you have an interesting marimo chemistry notebook, please submit a PR.

- smarts_view.py - a simple interactive viewer for SMILES and SMARTS
- reos_view.py - a viewer for functional group filters
- chem_view.py - interactive viewing of a projected t-SNE space
- ml_view.py - interactive viewing of machine learning model performance




