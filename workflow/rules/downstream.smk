rule apply_climate_change:
    input:
        tropcyc=tc_global_path("{track_model}", "no_cc"),
    output:
        tc_global_path("{track_model}", "{climate_scenario}"),
    wildcard_constraints:
        climate_scenario=r"(rcp26|rcp45|rcp60|rcp85)_(2100|20\d\d)",
    conda:
        "climada_env_TC"
    resources:
        mem_mb_per_cpu=24000,
    log:
        "logs/apply_climate_change_{track_model}_{climate_scenario}.log",
    script:
        "../scripts/apply_climate_change.py"


def get_all_countries(haz, countries):
    if countries == "all":
        import country_converter as coco
        from climada.hazard import TropCyclone

        haz = TropCyclone.from_hdf5(haz)
        haz.centroids.set_region_id()
        cc = coco.CountryConverter()
        found = cc.convert(
            haz.centroids.gdf["region_id"].unique(), to="ISO3", src="ISOnumeric"
        )
        countries = [c for c in found if c != "not found"]
    return countries


checkpoint tcs_per_countries:
    input:
        global_haz=TC_GLOBAL,
    output:
        directory("tropcyc/{track_model}/{climate_scenario}/by_country/"),
    params:
        countries=lambda w, input: get_all_countries(
            input.global_haz, config["countries"]
        ),
    conda:
        "climada_env_TC"
    resources:
        mem_mb_per_cpu=24000,
    log:
        "logs/get_countries_tcs_{track_model}_{climate_scenario}.log",
    script:
        "../scripts/get_country_tc.py"


def countries_of(track_model, climate_scenario):
    d = checkpoints.tcs_per_countries.get(
        track_model=track_model, climate_scenario=climate_scenario
    ).output[0]
    pattern = os.path.join(
        d, f"tropcyc_{tc_res_str}arcsec_{{country}}_{climate_scenario}.hdf5"
    )
    return glob_wildcards(pattern).country


def use_slr(wildcards):
    return wildcards.ssp != "nossp"


def get_mb_for_country_surge(wildcards, attempt):
    return max((attempt - 1) * 64000, 24000)


rule surges_countries:
    input:
        tropcyc=TC_COUNTRY,
        slr=branch(
            use_slr,
            then=expand(
                "{slr_data}/{{ssp}}/total_{{ssp}}_medium_confidence_values.nc",
                slr_data=config["slr_data_path"],
            ),
        ),
        dem=config["dem_topo_path"],
    output:
        SURGE_COUNTRY,
    params:
        higher_res=config["higher_res_surge"],
    conda:
        "climada_env_TC"
    resources:
        mem_mb_per_cpu=get_mb_for_country_surge,
        runtime=240,
    log:
        "logs/generate_surges_{track_model}_{country}_{ssp}_{slr_year}slr_{climate_scenario}.log",
    benchmark:
        "benchmark/generate_surge_{track_model}_{country}_{ssp}_{slr_year}slr_{climate_scenario}.tsv"
    script:
        "../scripts/generate_TCs_surges_country.py"


# ---- targets ----
rule all_tcs:
    input:
        [tc_global_path(m, s) for m, s in model_scenarios()],


def all_tcs_countries_inputs(wildcards):
    return [
        expand(
            TC_COUNTRY, track_model=m, climate_scenario=s, country=countries_of(m, s)
        )
        for m, s in model_scenarios()
    ]


rule all_tcs_countries:
    input:
        all_tcs_countries_inputs,


def all_surges_inputs(wildcards):
    slr = [(ssp, yr) for ssp in config["slr_ssp"] for yr in config["slr_years"]]
    slr.append(("nossp", "no"))
    files = []
    for m, s in model_scenarios():
        for country in countries_of(m, "no_cc"):
            files += (
                expand(
                    SURGE_COUNTRY,
                    track_model=m,
                    climate_scenario=s,
                    country=country,
                    ssp=[a for a, _ in slr][:0] or None,
                )
                if False
                else [
                    SURGE_COUNTRY.format(
                        track_model=m,
                        climate_scenario=s,
                        country=country,
                        ssp=ssp,
                        slr_year=yr,
                    )
                    for ssp, yr in slr
                ]
            )
    return files


rule all_surges_countries:
    input:
        all_surges_inputs,
