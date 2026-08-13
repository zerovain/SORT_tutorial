SORT tutorials
==============

This site contains tutorials for `SORT
<https://github.com/zerovain/SORT>`_ (Spatial Orthogonal Regularized
Transcriptomic decomposition), a method for discovering gene-expression
programs in spatial transcriptomic atlases.

.. image:: _static/sort_overview.png
   :alt: Overview of the SORT framework
   :width: 100%

Tutorials
---------

.. toctree::
   :maxdepth: 1

   Tutorial 1: Simulated spatial atlas <tutorials/00_simulated_atlas>
   Tutorial 2: PDAC atlas <tutorials/01_pdac_atlas>
   Tutorial 3: Mouse embryogenesis atlas <tutorials/02_embryogenesis>

Getting started
---------------

Install the SORT package first:

.. code-block:: bash

   git clone https://github.com/zerovain/SORT.git
   cd SORT
   conda env create -f environment.yml
   conda activate sort-paper-release
   python -m pip install -e . --no-deps

Then clone the tutorial repository:

.. code-block:: bash

   git clone https://github.com/zerovain/SORT_tutorial.git

The simulated input is included. The PDAC and mouse embryogenesis datasets can
be obtained from the public repositories linked in their respective tutorials.

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

Individual tutorials show the settings used for their datasets. Signal-profile
orthogonality is enforced directly by the Stiefel constraint throughout both
optimization stages.

Output
------

``result.W`` contains nonnegative spatial activities and ``result.Q`` contains
the corresponding gene weights shared across samples. The requested number of
components includes one background component. The first component (Python
index 0, reported as GEP00 in the manuscript) captures broad background
expression; the remaining components are the signal GEPs.

Citation
--------

Citation information will be provided upon publication of the associated
manuscript.

License
-------

The SORT package and tutorial code are released under the MIT License.
