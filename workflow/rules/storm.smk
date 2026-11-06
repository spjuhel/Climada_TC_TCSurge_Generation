rule storm_all_basins_no_cc:
    """STORM is present-climate only, so its output is directly the `no_cc` hazard."""
    input:
        expand(
            "tropcyc/storm/no_cc/by_basin/tropcyc_{tc_res_str}arcsec_{basin}_{ens}.hdf5",
            basin=TRACK_MODELS["storm"]["basins"],
            ens=TRACK_MODELS["storm"]["ensemble_members"],
            tc_res_str=tc_res_str,
        ),
    output:
        tc_global_path("storm", "no_cc"),
    conda:
        "../envs/climada_env.yaml"
    resources:
        mem_mb_per_cpu=24000,
    log:
        "logs/concatenate_storm_no_cc.log",
    script:
        "../scripts/concatenate_all_basins_tcs.py"


wildcard_constraints:
    basin=r"[A-Z]+",
    ens=r"\d+",


def storm_tracks_file(wildcards):
    """Path of the STORM track file for one basin and ensemble member."""
    return os.path.join(
        TRACK_MODELS["storm"]["tracks_dir"],
        f"STORM_DATA_IBTRACS_{wildcards.basin}_1000_YEARS_{wildcards.ens}.txt",
    )


rule storm_tcs_basin_ens:
    """Wind fields from STORM tracks for one basin and ensemble member"""
    input:
        tracks=storm_tracks_file,
        centroids="global_{tc_res_str}arcsec_centroids.hdf5",
    output:
        "tropcyc/storm/no_cc/by_basin/tropcyc_{tc_res_str}arcsec_{basin}_{ens}.hdf5",
    params:
        wind_model=TRACK_MODELS["storm"]["wind_model"],
        time_step_h=1.0,
        max_memory_gb=14,
    conda:
        "../envs/climada_env.yaml"
    resources:
        mem_mb_per_cpu=get_mb_for_split,
        runtime=480,
    log:
        "logs/storm_tcs_{tc_res_str}_{basin}_{ens}.log",
    benchmark:
        "benchmark/storm_tcs_{tc_res_str}_{basin}_{ens}.tsv"
    script:
        "../scripts/generate_TCs_storm_basin_ens.py"

rule storm_tcs_basin_ens_example:
    input:
        "tropcyc/storm/no_cc/by_basin/tropcyc_150arcsec_NI_0.hdf5",

rule storm_tcs_basin:
    """Wind fields from STORM tracks for one basin"""
    input:
        members = expand("tropcyc/storm/no_cc/by_basin/tropcyc_{{tc_res_str}}arcsec_{{basin}}_{ens}.hdf5", ens=TRACK_MODELS["storm"]["ensemble_members"])
    output:
        "tropcyc/storm/no_cc/by_basin/tropcyc_{tc_res_str}arcsec_{basin}.hdf5",
    conda:
        "../envs/climada_env.yaml"
    resources:
        mem_mb_per_cpu=get_mb_for_split,
    log:
        "logs/storm_tcs_{tc_res_str}_{basin}.log",
    benchmark:
        "benchmark/storm_tcs_{tc_res_str}_{basin}.tsv"
    script:
        "../scripts/generate_TCs_storm_basin_ens.py"

rule storm_tcs_basin_example:
    input:
        "tropcyc/storm/no_cc/by_basin/tropcyc_150arcsec_NI.hdf5",
