def tc_global_path(track_model, climate_scenario):
    """Standard path of the global hazard file (the model-agnostic interface)."""
    return (
        f"tropcyc/{track_model}/{climate_scenario}/"
        f"tropcyc_{tc_res_str}arcsec_all_basins_{climate_scenario}.hdf5"
    )


def model_scenarios():
    """All (track_model, climate_scenario) pairs requested in the config."""
    return [
        (m, s) for m, opts in TRACK_MODELS.items() for s in opts["climate_scenarios"]
    ]


# Pattern strings with wildcards, shared by downstream rules
TC_GLOBAL = tc_global_path("{track_model}", "{climate_scenario}")
TC_COUNTRY = (
    "tropcyc/{track_model}/{climate_scenario}/by_country/"
    f"tropcyc_{tc_res_str}arcsec_{{country}}_{{climate_scenario}}.hdf5"
)
SURGE_COUNTRY = (
    "surge/{track_model}/{climate_scenario}/by_country/"
    f"surge_{surge_res_str}arcsec_{{country}}_{{ssp}}_{{slr_year}}slr_{{climate_scenario}}.hdf5"
)

rule create_centroids:
    """Creates global centroids grid at `centroid_resolution` resolution"""
    output:
        "global_{tc_res_str}arcsec_centroids.hdf5",
    params:
        centroids_resolution=config["centroids_resolution"]/3600,
    log:
        "logs/create_centroids_{tc_res_str}.log",
    resources:
        mem_mb_per_cpu=12000,
    conda:
        "../envs/climada_env.yaml"
    script:
        "../scripts/create_centroids.py"
