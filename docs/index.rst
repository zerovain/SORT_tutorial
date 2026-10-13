SORT tutorials
==============

This site contains tutorials for `SORT
<https://github.com/zerovain/SORT>`_ (Spatial Orthogonal Regularized
Transcriptomic decomposition), a method for discovering gene-expression
programs in spatial transcriptomic atlases.

.. image:: _static/sort_overview.png
   :alt: Overview of the SORT framework
   :width: 100%

Getting started
---------------

Start with a recent Conda or Miniforge installation, then create the SORT
environment using the `installation guide
<https://github.com/zerovain/SORT#installation>`_:

The environment files supply Python 3.8.20 and the numerical core: NumPy, SciPy
and OpenBLAS. GPU acceleration is optional.

.. code-block:: bash

   git clone https://github.com/zerovain/SORT.git
   cd SORT
   conda env create --solver libmamba -f environment.yml

The default environment is tested on Linux x86_64. On Apple-silicon Macs,
replace ``environment.yml`` with ``environment-macos.yml`` in the command above.

Activate the environment:

.. code-block:: bash

   conda activate sort
   export PYTHONNOUSERSITE=1

**Apple silicon only:** build ``scikit-misc`` 0.2.0 before installing SORT.
With Apple's Command Line Tools and Homebrew installed:

.. code-block:: bash

   brew install gcc pkgconf
   export CC=/usr/bin/clang CXX=/usr/bin/clang++ FC=gfortran
   export SDKROOT="$(xcrun --show-sdk-path)"
   export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
   python -m pip install --no-cache-dir --no-deps --no-binary=scikit-misc scikit-misc==0.2.0

**Optional GPU acceleration.** Skip this step for CPU execution.
If you have an NVIDIA GPU, install CUDA-enabled
PyTorch and CuPy first. For a CUDA 12.1-compatible driver:

.. code-block:: bash

   python -m pip install torch==2.4.1 --index-url https://download.pytorch.org/whl/cu121 -c constraints-reference.txt
   python -m pip install cupy-cuda12x==12.3.0 -c constraints-reference.txt

See the `official installation guide <https://pytorch.org/get-started/locally/>`_
for platform and driver requirements. Keep PyTorch 2.4.1 for this environment.

**Install SORT.** Use the same command with or without GPU acceleration:

.. code-block:: bash

   python -m pip install -e '.[tutorial]' -c constraints-reference.txt
   python -m pip check

Conda installs the numerical dependencies. Create the environment before installing
SORT with pip to retain these numerical packages. The constraints fix the
remaining scientific dependencies, including Scanpy 1.9.8, Squidpy 1.2.3 and
``scikit-misc`` 0.2.0.

Use ``device="auto"`` when fitting to select GPU acceleration when available
or CPU execution otherwise.

Clone the tutorial repository and open the notebooks in JupyterLab:

.. code-block:: bash

   git clone https://github.com/zerovain/SORT_tutorial.git
   cd SORT_tutorial/notebooks
   jupyter lab

Choose the Python kernel from your SORT environment and run the cells in order.
Start with ``00_simulated_atlas.ipynb`` to fit SORT and explore its outputs
using the included data.

The simulated input is included. The PDAC and mouse embryogenesis datasets can
be obtained from the public repositories linked in their respective tutorials.
Each real-data tutorial describes the required AnnData input fields. The
simulation tutorial includes PCA/NMF comparisons, spatial maps and output exports.

Input data
----------

Use an AnnData object with locations as rows and genes as columns. Store gene
symbols in ``var_names``, spatial coordinates in ``obsm['spatial']`` and section
labels in ``obs``. For multiple sections, align genes to a common column order
and keep coordinates and metadata aligned with the expression rows. Coordinates
can remain in each section's own coordinate system.

The tutorials start from counts in ``adata.X`` and apply their stated
normalization and transformation before fitting. ``sort.fit`` takes the resulting
expression matrix. Build multi-section graphs with ``build_per_sample_graph``
using the section-label column, as illustrated in the real-data tutorials.

Core parameters
---------------

The notation below follows the SORT objective in the manuscript. The code
names are the arguments used by ``sort.fit`` or by the graph-building helper.

.. list-table:: Manuscript notation and Python arguments
   :header-rows: 1
   :widths: 18 25 57

   * - Manuscript
     - Python argument
     - Role
   * - :math:`K`
     - ``n_components``
     - Total factorization rank, including one background component and
       :math:`K-1` signal GEPs. Rank is dataset-specific; the manuscript uses
       the PCA scree-curve elbow as a guide.
   * - graph :math:`k`
     - ``n_neighbors`` in ``build_spatial_graph`` or
       ``build_per_sample_graph``
     - Number of spatial neighbors. Graphs for multi-section data are built
       within each section, without edges between sections.
   * - :math:`\lambda_{\mathrm{TV}}`
     - ``alpha``
     - Strength of graph total variation on ``W``. Larger values favor
       spatially coherent activities while retaining tissue boundaries.
   * - :math:`\lambda_{1,W}`
     - ``lambda_l1_W``
     - Base sparsity weight for spatial activities ``W``.
   * - :math:`\lambda_{1,Q}`
     - ``lambda_l1_Q``
     - Base sparsity weight for signal profiles ``Q``.
   * - :math:`\lambda_{-,Q}`
     - ``lambda_neg``
     - Penalizes negative values in the signal profiles.

The high-level defaults shared by most biological analyses are:

.. code-block:: python

   result = sort.fit(
       adata,
       n_components=25,         # choose for the dataset
       alpha=0.3,
       lambda_l1_W=0.3,
       lambda_l1_Q=300,
       lambda_neg=1.0,
   )

Individual tutorials show the settings used for their datasets. Signal profiles
remain orthogonal throughout both optimization stages.

Output
------

Rows of ``result.W`` correspond to locations, and rows of ``result.Q``
correspond to genes. Their columns identify the same components: ``W`` gives
nonnegative spatial activities and ``Q`` gives gene weights shared across samples.
The requested number of
components includes one background component. The first component (Python
index 0, reported as GEP00 in the manuscript) captures broad background
expression; the remaining components are the signal GEPs.

Inspect the genes with the largest positive weights to begin interpreting a
signal profile, then examine its spatial activity. For example:

.. code-block:: python

   import pandas as pd

   k = 1  # first signal component
   pd.Series(result.Q[:, k], index=result.gene_names).loc[
       lambda weights: weights > 0
   ].nlargest(10)

``result.save("outputs/my_atlas")`` exports the matrices, gene/location names
and metadata. Full dense reconstruction is opt-in; small slices can be computed
with ``result.reconstruct(rows=..., genes=...)``. When HVGs are selected, ``Q``
retains the input gene axis and ``result.metadata['fitted_gene_mask']``
identifies modeled genes.

Tutorials
---------

The example outputs shown in these tutorials were generated on Linux using
16 CPU threads; the real-data fits also used an NVIDIA A30 GPU.

.. toctree::
   :maxdepth: 1

   Tutorial 1: Simulated spatial atlas <tutorials/00_simulated_atlas>
   Tutorial 2: PDAC atlas <tutorials/01_pdac_atlas>
   Tutorial 3: Mouse embryogenesis atlas <tutorials/02_embryogenesis>

Citation
--------

Citation information will be provided upon publication of the associated
manuscript.

License
-------

The SORT package and tutorial code are released under the MIT License.
