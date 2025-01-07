# Tutorials


::::{grid} 1
:gutter: 1

:::{grid-item-card} Running on New Data {octicon}`database;1em;`
:link: data
:link-type: doc

Learn how to structure your input data, create configuration files, and run the pipeline on your own single-cell RNA sequencing datasets.
:::
:::{grid-item-card} Adding Methods {octicon}`code;1em;`
:link: method
:link-type: doc

Understand how to integrate your own transcriptional state discovery method into the CanSig benchmark framework for comparative analysis.
:::
:::{grid-item-card} Reproducing Results {octicon}`graph;1em;`
:link: reproducing
:link-type: doc

Step-by-step instructions for reproducing the benchmark results from our paper, including working with the curated cancer datasets.
:::
::::

## Overview of the snakemake pipeline

The diagram above shows the directed acyclic graph (DAG) of the Snakemake workflow, illustrating how different pipeline components connect and depend on each other.

![image](../assets/imgs/dag.pdf)



```{toctree}
:maxdepth: 2
:hidden: true
:titlesonly: true

data
method
reproducing
```