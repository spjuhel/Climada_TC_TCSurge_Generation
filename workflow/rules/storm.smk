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
        f"STORM_DATA_IBTRACS_{wildcards.genesis_basin}_1000_YEARS_{wildcards.ens}.txt",
    )

rule storm_tracks_basin:
    """Generate file with all tracks originating in `genesis_basin` for year `tracks_year`"""
    input:
        tracks=storm_tracks_file,
    output:
        "tracks/storm/{genesis_basin}/member_{ens}/STORM_{genesis_basin}_{ens}.hdf5",
    params:
        timestep=1,
    conda:
        "../envs/climada_env.yaml"
    resources:
        mem_mb_per_cpu=8000,
        runtime=240,
    log:
        "logs/find_tracks_storm_{genesis_basin}_{ens}.log",
    script:
        "../scripts/find_tracks_storm_basin.py"

checkpoint storm_split_basin_tracks:
    """Split tracks of a particular basin into files with at most `max_tracks_per_split` tracks"""
    input:
        tracks="tracks/storm/{genesis_basin}/member_{ens}/STORM_{genesis_basin}_{ens}.hdf5",
    output:
        directory("tracks/storm/{genesis_basin}/member_{ens}/splits/"),
    params:
        max_tracks=config["max_tracks_per_split"],
        output_name_format="STORM_{genesis_basin}_{ens}_split_{}.hdf5"
    conda:
        "../envs/climada_env.yaml"
    resources:
        mem_mb_per_cpu=8000,
        runtime=240,
    log:
        "logs/split_storm_tracks_{genesis_basin}_{ens}.log",
    script:
        "../scripts/split_storm_synth_tracks.py"


rule storm_tcs_basin_split_ens:
    """Wind fields from STORM tracks for one basin and ensemble member"""
    input:
        out_directory="tracks/storm/{genesis_basin}/member_{ens}/splits/",
        tracks="tracks/storm/{genesis_basin}/member_{ens}/splits/STORM_{genesis_basin}_{ens}_split_{split}.hdf5",
        centroids="global_{tc_res_str}arcsec_centroids.hdf5",
    output:
        "tropcyc/storm/no_cc/{genesis_basin}/member_{ens}/splits/STORM_{tc_res_str}arcsec_{genesis_basin}_{ens}_no_cc_split_{split}.hdf5",
    params:
        wind_model=TRACK_MODELS["storm"]["wind_model"],
        time_step_h=config["track_model"]["storm"]["time_step_h"],
        max_memory_gb=14,
        buffer=5,
    conda:
        "../envs/climada_env.yaml"
    resources:
        mem_mb_per_cpu=4000,
        runtime=480,
    log:
        "logs/storm_tcs_{tc_res_str}_{genesis_basin}_{ens}_{split}.log",
    benchmark:
        "benchmark/storm_tcs_{tc_res_str}_{genesis_basin}_{ens}_{split}.tsv"
    script:
        "../scripts/generate_TCs_storm_basin_ens.py"

def aggregate_storm_tcs(wildcards):
    checkpoint_output = checkpoints.storm_split_basin_tracks.get(**wildcards).output[0]
    return expand(
        "tropcyc/storm/no_cc/{{genesis_basin}}/member_{{ens}}/splits/STORM_{{tc_res_str}}arcsec_{{genesis_basin}}_{{ens}}_no_cc_split_{split}.hdf5",
        split=glob_wildcards(
            os.path.join(
                checkpoint_output,
                "STORM_{genesis_basin}_{ens}_split_{split}.hdf5"
            )
        ).split
    )

rule gather_split_storm_tcs:
    input:
        tracks=aggregate_storm_tcs,
    output:
        "tropcyc/storm/no_cc/{genesis_basin}/member_{ens}/STORM_{tc_res_str}arcsec_{genesis_basin}_{ens}_no_cc.hdf5",
    conda:
        "../envs/climada_env.yaml"
    resources:
        mem_mb_per_cpu=8000,
        runtime=240,
    log:
        "logs/gather_split_storm_TCs_{genesis_basin}_{ens}_{tc_res_str}arcsec.log",
    script:
        "../scripts/concatenate_all_splits_storm_tcs.py"


rule storm_tcs_basin_ens_example:
    input:
        "tropcyc/storm/no_cc/NI/member_0/STORM_600arcsec_NI_0_no_cc.hdf5",

rule storm_tcs_basin:
    """Wind fields from STORM tracks for one basin"""
    input:
        tracks = expand("tropcyc/storm/no_cc/{{genesis_basin}}/member_{ens}/STORM_{{tc_res_str}}arcsec_{{genesis_basin}}_{ens}_no_cc.hdf5", ens=config["track_model"]["storm"]["ensemble_members"])
    output:
        "tropcyc/storm/no_cc/{genesis_basin}/STORM_{tc_res_str}arcsec_{genesis_basin}_no_cc.hdf5",
    conda:
        "../envs/climada_env.yaml"
    resources:
        mem_mb_per_cpu=get_mb_for_split,
    log:
        "logs/storm_tcs_{tc_res_str}_{genesis_basin}.log",
    benchmark:
        "benchmark/storm_tcs_{tc_res_str}_{genesis_basin}.tsv"
    script:
        "../scripts/concatenate_all_splits_storm_tcs.py"

rule storm_tcs_basin_example:
    input:
        "tropcyc/storm/no_cc/NI/STORM_600arcsec_NI_no_cc.hdf5",

rule storm_tcs:
    input:
        expand("tropcyc/storm/no_cc/{genesis_basin}/STORM_{tc_res_str}arcsec_{genesis_basin}_no_cc.hdf5",
               genesis_basin=config["track_model"]["storm"]["basins"],
               tc_res_str=config["centroids_resolution"]
               )
