SORT tutorials
==============

This site contains tutorials for `SORT
<https://github.com/zerovain/SORT>`_ (Spatial Orthogonal Regularized
Transcriptomic decomposition), a method for discovering gene-expression
programs in spatial transcriptomic atlases.

Tutorials
---------

.. toctree::
   :maxdepth: 1

   Tutorial 1: Simulated spatial atlas <tutorials/00_simulated_atlas>
   Tutorial 2: PDAC atlas <tutorials/01_pdac_atlas>
   Tutorial 3: Mouse embryogenesis atlas <tutorials/02_embryogenesis>

.. toctree::
   :hidden:

   TUTORIALS
   DATA_AVAILABILITY

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

The simulated input is included. The PDAC and mouse embryogenesis tutorials
provide preparation scripts for data obtained from the original public
repositories; processed third-party datasets are not redistributed.

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
