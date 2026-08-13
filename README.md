# SORT tutorials

Tutorials for [SORT](https://github.com/zerovain/SORT), a spatial
transcriptomics method for learning spatially coherent gene-expression
programs.

## Installation

Clone both repositories into the same parent directory and create the tested
environment:

```bash
git clone https://github.com/zerovain/SORT.git
git clone https://github.com/zerovain/SORT_tutorial.git
conda env create -f SORT/environment.yml
conda activate sort-paper-release
python -m pip install -e 'SORT[docs]'
```

Open the notebooks with Jupyter:

```bash
jupyter lab SORT_tutorial/notebooks
```

## Tutorials

1. [Simulated spatial atlas](notebooks/00_simulated_atlas.ipynb): complete
   seed-03 simulation with 10,368 locations, 5,000 modeled genes and eight
   sections.
2. [Complete PDAC atlas](notebooks/01_pdac_atlas.ipynb): complete 108-section
   GSE282302 PDAC atlas.
3. [Mouse embryogenesis atlas](notebooks/02_embryogenesis.ipynb): complete
   eight-stage MOSTA embryogenesis series.

Only the author-generated simulation H5AD is distributed. The biological
tutorials link to the original repositories and use the assembly scripts under
`prepare/`; processed third-party H5AD files stay local. Normalization,
`log1p`, highly variable gene selection and graph construction are shown in the
notebooks.

The paired Python sources can regenerate the notebooks or export static pages:

```bash
python build_docs.py
python build_docs.py --format html
```

Full paper reproduction, including datasets not used as public tutorials,
remains in the companion manuscript-analysis archive.

## License

The tutorial code is released under the [MIT License](LICENSE).
