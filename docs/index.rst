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
   * - :math:`\varepsilon_{\mathrm{TV}}`
     - ``tv_epsilon``
     - Smooths the absolute differences used by graph TV for numerical
       stability.
   * - Two optimization stages
     - ``stage1_epochs``, ``stage2_epochs`` and ``adam_steps``
     - Sets the Stage 1 and Stage 2 epoch counts and the number of Riemannian
       Adam profile updates per Stage 2 epoch.
   * - TV update schedule
     - ``tv_stage`` and ``tv_update_freq``
     - Selects the stage in which graph TV is active and how often its IRLS
       weights are refreshed.
   * - random seed
     - ``random_state``
     - Controls reproducible initialization.

With ``l1_weight_strategy="adaptive"``, ``lambda_l1_W`` and
``lambda_l1_Q`` are base values used to calculate component-specific sparsity
weights after initialization; those weights are then held fixed during fitting.
The high-level defaults shared by most biological analyses are:

.. code-block:: python

   result = sort.fit(
       adata,
       n_components=25,         # choose for the dataset
       alpha=0.3,
       lambda_l1_W=0.3,
       lambda_l1_Q=300,
       lambda_neg=1.0,
       use_tv=True,
       tv_epsilon=0.01,
       tv_update_freq=5,
       tv_stage="stage2",
       stage1_epochs=50,
       stage2_epochs=100,
       adam_steps=40,
       random_state=42,
   )

Individual tutorials show the settings used for their datasets. ``beta`` is a
deprecated argument retained only for compatibility with older scripts; it is
not part of the current loss and should not be tuned. Signal-profile
orthogonality is enforced directly by the Stiefel constraint throughout both
optimization stages. The internal Adam coefficients ``beta1`` and ``beta2``
are unrelated to this deprecated argument.

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
