from itertools import product

from config_parser import ParsedConfig

import os
import sys

# Add project root to path
workflow_dir = workflow.basedir
sys.path.insert(0, workflow_dir)
os.environ['PYTHONPATH'] = f"{workflow_dir}/scripts:{os.environ.get('PYTHONPATH', '')}"
print(os.environ['PYTHONPATH'] )
cfg = ParsedConfig(config)


include: "rules/preprocessing/Snakefile"
include: "rules/integration/Snakefile"
include: "rules/signatures/Snakefile"
include: "rules/evaluate_signatures/Snakefile"

localrules: all

wildcard_constraints:
    type=r"scrna|bulk|overlap",
    method=r"|".join(cfg.get_methods()),
    signature=r"|".join(cfg.signatures),
    n_samples=r"\d+|all",
    random_seed = r"\d+",
    n_hvg = r"\d+|all",
    batch_key=r"True|False"

rule all:
    input: expand(expand(rules.aggregate_signatures.output, zip, **cfg.get_scenarios_types(), allow_missing=True), n_samples=cfg.n_samples, n_hvg=cfg.n_hvgs, batch_key=cfg.use_batch_key)


