import pathlib
from collections import defaultdict
from itertools import product


class ParsedConfig:
    def __init__(self, config):
        self.methods = config["methods"]
        self.ROOT = pathlib.Path(config["root"]).resolve()
        self.signatures = config["signatures"]
        self.scenarios = config["scenarios"]
        self.validate_scenarios()

    def validate_scenarios(self):
        for scenario in self.scenarios.values():
            for signature in scenario["signatures"]:
                if signature not in self.signatures:
                    raise ValueError(f"Unknown signature {signature} in {scenario}")

    def get_scenarios(self):
        return list(self.scenarios.keys())

    def get_preprocessing_dict(self, wildcards):
        if "preprocessing" in self.scenarios[wildcards.scenario]:
            return self.scenarios[wildcards.scenario]["preprocessing"]
        return {}

    def get_excluded_samples(self, wildcards):
        config = self.get_preprocessing_dict(wildcards)
        return config.get("excluded_samples", None)

    def get_min_genes(self, wildcards):
        config = self.get_preprocessing_dict(wildcards)
        return config.get("min_genes", 1000)

    def get_max_pct_mt(self, wildcards):
        config = self.get_preprocessing_dict(wildcards)
        return config.get("max_pct_mt", 20)

    def get_scenarios_types(self):
        scenarios_types = defaultdict(list)
        for scenario in self.get_scenarios():
            scenarios_types["scenario"].append(scenario)
            scenarios_types["type"].append("scrna")
            scenarios_types["scenario"].append(scenario)
            scenarios_types["type"].append("overlap")
            if "bulk_path" in self.scenarios[scenario].keys():
                scenarios_types["scenario"].append(scenario)
                scenarios_types["type"].append("bulk")
        return scenarios_types

    def get_signatures_for_scenario(self, wildcards):
        return self.scenarios[wildcards.scenario]["signatures"]

    def get_data_path(self, wildscards):
        if wildscards.scenario in self.get_scenarios():
            return self.scenarios[wildscards.scenario]["data_path"]
    def get_methods(self):
        return list(self.methods.keys())
    
    def get_cansig_metasig_method(self, wildscards):
        if not wildscards.method.startswith("cansig"):
            raise ValueError(f"This method is only for the use with CanSig but was used with {wildscards.method}")
        return self.methods[wildscards.method]["metasig_method"]

    def get_n_cluster(self, wildcards):
        return self.signatures[wildcards.signature]["n_cluster"]

    def get_bulk_path(self, wildcards):
        return self.scenarios[wildcards.scenario]["bulk_path"]

    def get_scrna_path(self, wildcards):
        return self.signatures[wildcards.signature]["scoring_path"]

    def get_scoring_scenario(self, wildcards):
        return self.signatures[wildcards.signature]["scoring_scenario"]

    def get_options(self, wildcards):
        config = self.methods[wildcards.method]
        if config:
            return " ".join([f"--{key} {value}" for key, value in config.items()])
        return ""

    def get_annotation_path(self, wildcards):
        return self.signatures[wildcards.signature]["annotation_path"]

    def get_n_samples(self, wildcards) -> list:
        scenario = self.scenarios[wildcards.scenario]
        if "n_samples" in scenario:
            return scenario["n_samples"]
        else:
            return ["all"]

    def get_random_seeds(self, wildcards) -> list:
        scenario = self.scenarios[wildcards.scenario]
        if "random_seeds" in scenario:
            return list(range(scenario["random_seeds"]))
        else:
            return [0]