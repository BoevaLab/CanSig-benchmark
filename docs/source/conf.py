# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
from pathlib import Path


root_dir = Path(__file__).parents[2].absolute()
print(root_dir)
sys.path.insert(0, str(root_dir))

project = 'CanSig'
copyright = '2024, Florian Barkmann, Josephine Yates'
author = 'Florian Barkmann, Josephine Yates'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['myst_parser', 
              "sphinx_design",
            'sphinx.ext.duration',
            'sphinx.ext.doctest',
            'sphinx.ext.autodoc',
            'sphinx.ext.autosummary',
            'sphinx.ext.napoleon',
            'sphinx.ext.viewcode' 
              ]
myst_enable_extensions = ["colon_fence"]


templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
html_logo = 'assets/imgs/logo.png'

pygments_style = "default"

import os
import sys

for path in ["preprocessing", "metasig", "metrics"]:
  sys.path.insert(0, os.path.abspath(f'scripts'))