def tc_global_path(track_model, climate_scenario):
    """Standard path of the global hazard file (the model-agnostic interface)."""
    return (
        f"tropcyc/{track_model}/{climate_scenario}/"
        f"tropcyc_{tc_res_str}arcsec_all_basins_{climate_scenario}.hdf5"
    )

def model_scenarios():
    """All (track_model, climate_scenario) pairs requested in the config."""
    return [
        (m, s)
        for m, opts in TRACK_MODELS.items()
        for s in opts["climate_scenarios"]
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
