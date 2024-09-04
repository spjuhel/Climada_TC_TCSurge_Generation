import os
import sys
import logging, traceback
from pathlib import Path

import numpy as np
import re

from climada.hazard import TCTracks, TropCyclone, Centroids

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )



# Install exception handler
sys.excepthook = handle_exception


climate_scenarios = snakemake.config["climate_scenarios"]
climate_sce_re = re.compile(r"no_cc|rcp(\d)(\d)_(2100|20\d\d)")
tropcyc = snakemake.input.tropcyc
tc_res_str = snakemake.params.tc_res_str
nsynth = snakemake.params.nsynth
start_period = snakemake.params.start_period
end_period = snakemake.params.end_period


logger.info(f"Reading TC events from {tropcyc}")
if os.stat(tropcyc).st_size == 0:
    logger.info(
        f"File is empty, which probably means there is no TCs data. Ignoring."
    )
    for climate_scenario in climate_scenarios:
        out = f"tropcyc/{climate_scenario}/tropcyc_{tc_res_str}arcsec_{nsynth}synth_all_basins_{start_period}_to_{end_period}_{climate_scenario}.hdf5"
        logger.info(f"Writing to {out}")
        Path(out).touch()
else:
    tc = TropCyclone.from_hdf5(tropcyc)
    for climate_scenario in climate_scenarios:
        scenario = climate_sce_re.match(climate_scenario)
        if not scenario:
            raise ValueError(f"Not a valid climate scenario: {climate_scenario}")
        if climate_scenario != "no_cc":
            rcp_arg = f"{scenario.group(1)}.{scenario.group(2)}"
            cc_ref_year = int(scenario.group(3))
            logger.info(f"Applying climate change (rcp{rcp_arg} - {cc_ref_year})")
            tc_clim = tc.apply_climate_scenario_knu(target_year=int(cc_ref_year), scenario=rcp_arg)
            out = f"tropcyc/{climate_scenario}/tropcyc_{tc_res_str}arcsec_{nsynth}synth_all_basins_{start_period}_to_{end_period}_{climate_scenario}.hdf5"
            logger.info(f"Writing to {out}")
            tc_clim.write_hdf5(out)
