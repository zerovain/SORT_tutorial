# SORT tutorials

This repository contains tutorials for
[SORT](https://github.com/zerovain/SORT) (Spatial Orthogonal Regularized
Transcriptomic decomposition), a computational method for discovering
gene-expression programs in spatial transcriptomic atlases.

![Overview of the SORT framework](docs/_static/sort_overview.png)

## Contents

- Tutorial 1: [Simulated spatial atlas](notebooks/00_simulated_atlas.ipynb)
- Tutorial 2: [PDAC atlas](notebooks/01_pdac_atlas.ipynb)
- Tutorial 3: [Mouse embryogenesis atlas](notebooks/02_embryogenesis.ipynb)

The simulated dataset is included in this repository. The PDAC and mouse
embryogenesis datasets can be obtained from the public repositories linked in
their respective tutorials.

## Getting started

Install SORT by following the instructions in the
[SORT repository](https://github.com/zerovain/SORT#installation), then clone
this repository:

```bash
git clone https://github.com/zerovain/SORT_tutorial.git
```

The tutorials can be viewed directly on GitHub or opened locally with Jupyter:

```bash
jupyter lab SORT_tutorial/notebooks
```

## Documentation

Online documentation will be made available through Read the Docs. Until then,
the notebooks in this repository contain the complete tutorials.

## License

This project is licensed under the [MIT License](LICENSE).
