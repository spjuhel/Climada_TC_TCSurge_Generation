rule tracks_basin_year:
    """Generate file with all tracks originating in `genesis_basin` for year `tracks_year`"""
    output:
        "tracks/{genesis_basin}/{tracks_year}/IBTracs_{genesis_basin}_{tracks_year}.hdf5",
    params:
        timestep=1,
    conda:
        "climada_env_TC"
    resources:
        mem_mb_per_cpu=8000,
        runtime=240,
    log:
        "logs/find_tracks_{genesis_basin}_{tracks_year}.log",
    script:
        "../scripts/find_tracks_basin_year.py"


def synth_track_memb(wildcards, attempt):
    return 1000 * (config["nsynth"] // 2) * attempt


rule synth_tracks_basin_year:
    """Generate file with `nsynth` synthetic tracks from historical ones in `genesis_basin` for year `tracks_year`"""
    input:
        "tracks/{genesis_basin}/{tracks_year}/IBTracs_{genesis_basin}_{tracks_year}.hdf5",
    output:
        f"tracks/{{genesis_basin}}/{{tracks_year}}/IBTracs_{nsynth}synth_{{genesis_basin}}_{{tracks_year}}.hdf5",
    params:
        nsynth=nsynth,
    conda:
        "climada_env_TC"
    resources:
        mem_mb_per_cpu=synth_track_memb,
        runtime=240,
    log:
        "logs/generate_synths_tracks_{genesis_basin}_{tracks_year}.log",
    script:
        "../scripts/generate_synths_tracks_basin_year.py"


# Naming is not automagic in the script!!!
checkpoint split_basin_year_tracks:
    """Split synthetic tracks of a particular basin, year into files with at most `max_tracks_per_split` tracks"""
    input:
        tracks=f"tracks/{{genesis_basin}}/{{tracks_year}}/IBTracs_{nsynth}synth_{{genesis_basin}}_{{tracks_year}}.hdf5",
    output:
        directory("tracks/{genesis_basin}/{tracks_year}/synth_splits/"),
    params:
        max_tracks=config["max_tracks_per_split"],
    conda:
        "climada_env_TC"
    resources:
        mem_mb_per_cpu=8000,
        runtime=240,
    log:
        "logs/split_synth_tracks_{genesis_basin}_{tracks_year}.log",
    script:
        "../scripts/split_synth_tracks.py"


def get_mb_for_split(wildcards, attempt):
    return max((attempt - 1) * 64000, 20000)


rule tcs_basin_year_split:
    """Generate tropical cyclones from synthetic tracks (for a split file)"""
    input:
        tracks=f"tracks/{{genesis_basin}}/{{tracks_year}}/synth_splits/IBTracs_{nsynth}synth_{{genesis_basin}}_{{tracks_year}}_split_{{split}}.hdf5",
        centroids=f"global_{tc_res_str}arcsec_centroids.hdf5",
    output:
        f"tropcyc/{{genesis_basin}}/no_cc/{{tracks_year}}/tropcyc_{tc_res_str}arcsec_{nsynth}synth_{{genesis_basin}}_{{tracks_year}}_no_cc_split_{{split}}.hdf5",
    params:
        buf=5,
        max_memory_gb=14,
    conda:
        "climada_env_TC"
    resources:
        mem_mb_per_cpu=get_mb_for_split,
        runtime=480,
    log:
        "logs/generate_TCs_{genesis_basin}_{tracks_year}_split_{split}.log",
    benchmark:
        "benchmark/generate_TCs_{genesis_basin}_{tracks_year}_split_{split}.tsv"
    script:
        "../scripts/generate_TCs_basin_year.py"


def aggregate_tcs(wildcards):
    checkpoint_output = checkpoints.split_basin_year_tracks.get(**wildcards).output[0]
    return expand(
        "tropcyc/{{genesis_basin}}/no_cc/{{tracks_year}}/tropcyc_{tc_res_str}arcsec_{nsynth}synth_{{genesis_basin}}_{{tracks_year}}_no_cc_split_{split}.hdf5",
        split=glob_wildcards(
            os.path.join(
                checkpoint_output,
                f"IBTracs_{nsynth}synth_{{genesis_basin}}_{{tracks_year}}_split_{{split}}.hdf5",
            )
        ).split,
        nsynth=nsynth,
        tc_res_str=tc_res_str,
    )


rule gather_split_tcs:
    input:
        aggregate_tcs,
    output:
        f"tropcyc/{{genesis_basin}}/no_cc/{{tracks_year}}/tropcyc_{tc_res_str}arcsec_{nsynth}synth_{{genesis_basin}}_{{tracks_year}}_no_cc.hdf5",
    conda:
        "climada_env_TC"
    resources:
        mem_mb_per_cpu=16000,
        runtime=240,
    log:
        "logs/gather_split_TCs_{genesis_basin}_{tracks_year}_historical.log",
    script:
        "../scripts/concatenate_all_basins_tcs.py"


rule concatenate_all_years_tcs:
    input:
        expand(
            "tropcyc/{{genesis_basin}}/no_cc/{tracks_year}/tropcyc_{tc_res_str}arcsec_{nsynth}synth_{{genesis_basin}}_{tracks_year}_no_cc.hdf5",
            tracks_year=tracks_year_range,
            tc_res_str=tc_res_str,
            nsynth=nsynth,
        ),
    output:
        f"tropcyc/{{genesis_basin}}/no_cc/tropcyc_{tc_res_str}arcsec_{nsynth}synth_{{genesis_basin}}_{start_period}_to_{end_period}_no_cc.hdf5",
    resources:
        mem_mb_per_cpu=24000,
    conda:
        "climada_env_TC"
    log:
        "logs/concatenate_all_years_tcs_{genesis_basin}_historical.log",
    script:
        "../scripts/concatenate_all_years_tcs.py"


rule all_no_cc_tcs:
    input:
        f"tropcyc/no_cc/tropcyc_{tc_res_str}arcsec_{nsynth}synth_all_basins_{start_period}_to_{end_period}_no_cc.hdf5",


rule concatenate_all_basins_tcs:
    input:
        expand(
            "tropcyc/{genesis_basin}/no_cc/tropcyc_{tc_res_str}arcsec_{nsynth}synth_{genesis_basin}_{start_period}_to_{end_period}_no_cc.hdf5",
            genesis_basin=TRACK_MODELS["ibtracs"]["genesis_basins"],
            tc_res_str=tc_res_str,
            nsynth=nsynth,
            start_period=start_period,
            end_period=end_period,
        ),
    output:
        tc_global_path("ibtracs", "no_cc"),
    conda:
        "climada_env_TC"
    resources:
        mem_mb_per_cpu=24000,
    log:
        "logs/concatenate_all_basins_tcs_ibtracs_no_cc.log",
    script:
        "../scripts/concatenate_all_basins_tcs.py"
