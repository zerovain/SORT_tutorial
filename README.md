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

Install SORT with the common installation command in the
[SORT repository](https://github.com/zerovain/SORT#installation), including
the `tutorial` extra for JupyterLab and `constraints-reference.txt` to select
the tested scientific dependencies. Create the Python 3.8.20 environment
from the SORT repository using `environment.yml` (tested on Linux x86_64).
On Apple-silicon Macs, use `environment-macos.yml` and complete the
`scikit-misc` 0.2.0 build step in the installation guide before installing SORT.
If you have an NVIDIA GPU, complete the optional PyTorch/CuPy step before
installing SORT. Otherwise, skip that step. Use `device="auto"` when fitting to select
the available device. Then clone this repository:

```bash
git clone https://github.com/zerovain/SORT_tutorial.git
```

The tutorials can be viewed directly on GitHub or opened locally with Jupyter:

```bash
cd SORT_tutorial/notebooks
jupyter lab
```

Use the Python kernel from the environment where SORT was installed. The
simulation notebook includes the input data, PCA/NMF comparisons and output
exports. Each real-data notebook provides the public data source and the
required H5AD input fields.

The example outputs shown in these tutorials were generated on Linux using
16 CPU threads; the real-data fits also used an NVIDIA A30 GPU.

## Input data

Use AnnData with locations as rows, genes as columns, gene symbols in
`var_names`, coordinates in `obsm['spatial']` and section labels in `obs`.
Align genes to a common column order across sections, keeping expression,
coordinates and metadata aligned by row. Coordinates can remain in each
section's own coordinate system.

The tutorials start from counts in `adata.X`, then normalize and transform
expression before passing it to `sort.fit`. Build graphs within each section
using `build_per_sample_graph` and the section-label column.

PDAC uses `obs['sample_id']` for section identifiers; embryogenesis uses
`obs['slice_id']` and `obs['annotation']` for stage and anatomical labels,
with one section per stage in this example. Select components in your fit
using their gene profiles and spatial patterns; the indices shown refer to
the example outputs.

## Documentation

The full documentation is available at
[sort-tutorial.readthedocs.io](https://sort-tutorial.readthedocs.io/).

## License

This project is licensed under the [MIT License](LICENSE).
