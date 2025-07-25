from config_parser import ParsedConfig

cfg = ParsedConfig(config)

workflow_dir = workflow.basedir
sys.path.insert(0, workflow_dir)
os.environ['PYTHONPATH'] = f"{workflow_dir}/scripts:{os.environ.get('PYTHONPATH', '')}"
print(os.environ['PYTHONPATH'] )

include: "rules/preprocessing/Snakefile"
include: "rules/integration/Snakefile"
include: "rules/evaluate_integration/Snakefile"

localrules: all

wildcard_constraints:
    type="scrna|bulk|overlap",
    method="|".join(cfg.get_methods()),
    n_samples=r"all|-?\d+\.\d+|-?\d+",  
    random_seed=r"\d+"

print("DEBUG INFO:")
print(f"Scenarios: {cfg.get_scenarios()}")
print(f"Methods: {cfg.get_methods()}")
print(f"n_hvgs: {cfg.n_hvgs}")
print(f"batch_key: {cfg.use_batch_key}")
print(f"ROOT: {cfg.ROOT}")


rule all:
    input:
        metrics=expand(
            rules.eval_latent.output.metrics, 
            scenario=cfg.get_scenarios(), 
            method=cfg.get_methods(), 
            random_seed=cfg.random_seed, 
            n_samples=cfg.n_samples, 
            n_hvg=cfg.n_hvgs, 
            batch_key=cfg.use_batch_key
        ),
        plotting=expand(
            rules.eval_latent.output.umap, 
            scenario=cfg.get_scenarios(),
            method=cfg.get_methods(),
            random_seed=0, 
            n_samples="all", 
            n_hvg=cfg.n_hvgs, 
            batch_key=cfg.use_batch_key
        ) if cfg.plotting else []