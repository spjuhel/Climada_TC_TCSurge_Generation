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


basin = snakemake.wildcards.genesis_basin
climate_scenarios = snakemake.config["climate_scenarios"]
climate_sce_re = re.compile(r"historical|rcp(\d)(\d)_(2100|20\d\d)")
tropcyc = snakemake.input.tropcyc

logger.info(f"Applying climate change scenario to genesis basin {basin}")

logger.info(f"Reading TC events from {tropcyc}")
tc = TropCyclone.from_hdf5(tropcyc)


for climate_scenario in climate_scenarios:
    scenario = climate_sce_re.match(climate_scenario)
    if not scenario:
        raise ValueError(f"Not a valid climate scenario: {climate_scenario}")
    if climate_scenario != "historical":
        rcp_arg = f"{scenario.group(1)}.{scenario.group(2)}"
        cc_ref_year = int(scenario.group(3))
        logger.info(f"Applying climate change (rcp{rcp_arg} - {cc_ref_year})")
        tc_clim = tc.apply_climate_scenario_knu(target_year=int(cc_ref_year), scenario=rcp_arg)
        out = f"tropcyc/{basin}/{climate_scenario}/TCs_{basin}_all_{climate_scenario}.hdf5"
        logger.info(f"Writing to {out}")
        tc_clim.write_hdf5(out)
