API Reference
=============

.. note:: 
    The documentation for R modules will be added in the future. Please reach out via github if you have any questions regading them.

The API consists of three main modules:

.. grid:: 1
   :gutter: 1

   .. card:: Preprocessing module |database|
       :link: preprocessing
       :link-type: doc

       This module handles the preprocessing and subsampling of the data.

   .. card:: Meta-signature module |code|
       :link: metasigs
       :link-type: doc
       
       This module generates gene signatures using different early and late integration methods.

   .. card:: Metrics module |graph| 
       :link: metrics
       :link-type: doc

       This module evaluates the similarities of the generated signatures with the provided signatures.

.. |database| replace:: :octicon:`database;1em`
.. |code| replace:: :octicon:`code;1em`
.. |graph| replace:: :octicon:`graph;1em`


.. toctree::
   :maxdepth: 2
   :hidden:

   preprocessing
   metasigs
   metrics
