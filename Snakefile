from itertools import product

from config_parser import ParsedConfig

cfg = ParsedConfig(config)


include: "scripts/preprocessing/Snakefile"
include: "scripts/metasigs/scalop/Snakefile"
include: "scripts/metasigs/early_integration/Snakefile"
include: "scripts/metrics/Snakefile"
include: "scripts/metasigs/genenmf/Snakefile"

localrules: all


wildcard_constraints:
    type="scrna|bulk|overlap",
    method="|".join(cfg.get_methods()),
    n_samples="\d+|all",
    random_seed = "\d+"

rule all:
    input: expand(rules.aggregate_signatures.output, zip, **cfg.get_scenarios_types())
