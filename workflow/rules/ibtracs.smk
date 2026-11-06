rule concatenate_all_basins_tcs:
    input:
        expand(
            "tropcyc/{genesis_basin}/no_cc/tropcyc_{tc_res_str}arcsec_{nsynth}synth_{genesis_basin}_{start_period}_to_{end_period}_no_cc.hdf5",
            genesis_basin=TRACK_MODELS["ibtracs"]["genesis_basins"],
            tc_res_str=tc_res_str, nsynth=nsynth,
            start_period=start_period, end_period=end_period,
        ),
    output:
        tc_global_path("ibtracs", "no_cc"),
    conda: "climada_env_TC"
    resources: mem_mb_per_cpu=24000
    log: "logs/concatenate_all_basins_tcs_ibtracs_no_cc.log"
    script: "../scripts/concatenate_all_basins_tcs.py"
