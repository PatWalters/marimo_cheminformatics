# marimo_cheminformatics
Practical Cheminformatics with Marimo

The notebooks in this repo show some examples of how to use marimo notebooks as a tool for cheminformatic analysis. Marimo is often referred to as a "better Jupyter notebook". I've found two primary advantages of marimo over Jupyter.
- Marimo is reactive, changes to a variable in one cell will be updated all other cells containing that variable.
- Marimo contains several native GUI components that can be integrated to easily generate interactive analyses and dashboards.

For more information on marimo, I'd recommend these tutorials and references.
- [Marimo's getting started docs](https://docs.marimo.io/getting_started/)
- [Marimo FAQ](https://docs.marimo.io/faq/)
- [Marimo article in Real Python](https://realpython.com/marimo-notebook/)

For those who prefer videos, these are excellent.
- [An overview of marimo](https://www.youtube.com/watch?v=3N6lInzq5MI)
- [Python notebooks: Marimo vs Jupyter](https://www.youtube.com/watch?v=tLyjRfkyfFg&t=265s)
- [Marimo Notebooks Intro | Charting Python's rise in popularity](https://www.youtube.com/watch?v=XoArtLKPJ2I)

### Installation
1. Clone this repo.
```git clone git@github.com:PatWalters/marimo_cheminformatics.git```
2. Create an environment and install the necessary libraries.
Not using [uv](https://docs.astral.sh/uv/) yet? Here's a [great tutorial](https://realpython.com/python-uv/) explaining why you should. 
```uv venv marimo_cinf --python 3.11
source marimo_conf/bin/activate
cd marimo_marimocheminformatics
uv pip install -e .
```
3. Install marimo
The marimo docs provide example of how to install marimo with [pip](https://docs.marimo.io/#__tabbed_1_1),
[uv](https://docs.marimo.io/#__tabbed_1_2), or [conda](https://docs.marimo.io/#__tabbed_1_3).
All you really have to do is this.    
```uv pip install marimo```

### Running the notebooks
The notebooks can be run with `marimo edit notebook.py`, where `notebook.py` is the placehold for the notebook you want to run.

### The current collection
This repo currently contains a few notebooks demonstrating how to use marimo for cheminformatics. I'll continue add notebooks. If you have an interesting marimo chemistry notebook, please submit a PR.

- smarts_view.py - a simple interactive viewer for SMILES and SMARTS
- reos_view.py - a viewer for functional group filters
- chem_view.py - interactive viewing of a projected t-SNE space
- ml_view.py - interactive viewing of machine learning model performance

### Changes and fixes are appreciated
If you come across and error or a typo, please file an issue or submit a PR.  In addtion, if you have improvements to these notebooks, please submit a PR. I'd love it if we could make this repo a place to share how marimo can facilitate cheminformatics analyses.







