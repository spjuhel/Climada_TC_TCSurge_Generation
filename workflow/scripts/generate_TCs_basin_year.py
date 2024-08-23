import os
import sys
import logging, traceback
from pathlib import Path
import re

import numpy as np

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
year = snakemake.wildcards.tracks_year
split = snakemake.wildcards.split
climate_scenarios = snakemake.config.climate_scenarios
climate_sce_re = re.compile(r"historical|rcp(\d)(\d)_(2100|20\d\d)")
logger.info(f"Computing TC events for genesis basin {basin} for {year}")

logger.info(f"Loading TC tracks from {snakemake.input.tracks}")

if os.stat(snakemake.input[0]).st_size == 0:
    logger.info(f"File is empty, which probably means there is no track data for this basin-year. Ignoring")
    Path(snakemake.output[0]).touch()
else:
    tracks = TCTracks.from_hdf5(snakemake.input.tracks)
    logger.info(f"There are {len(tracks.data)} tracks.")

    logger.info(f"Loading global centroids from {snakemake.input.centroids}")
    cent = Centroids.from_hdf5(snakemake.input.centroids)

    logger.info(f"Selecting centroids extent from tracks with buffer={snakemake.params.buf}")
    cent_tracks = cent.select(extent=tracks.get_extent(snakemake.params.buf))

    logger.info(f"Computing TC wind-fields")
    tracks.equal_timestep(0.1)
    tc = TropCyclone.from_tracks(tracks, centroids=cent_tracks, max_memory_gb=snakemake.params.max_memory_gb)
    freq_corr = 1 / snakemake.config["nsynth"]
    tc.frequency = np.ones(tc.event_id.size)*freq_corr
    out = f"tropcyc/{basin}/TCs_{basin}_{year}_historical_split_{split}.hdf5"
    logger.info(f"Writing to {out}")
    tc.write_hdf5(out)
    for climate_scenario in climate_scenarios:
        scenario = climate_sce_re.match(climate_scenario)
        if not scenario:
            raise ValueError(f"Not a valid climate scenario: {climate_scenario}")
        if climate_scenario != "historical":
            rcp_arg = f"{scenario.group(1)}.{scenario.group(2)}"
            cc_ref_year = int(scenario.group(3))
            logger.info(f"Applying climate change (rcp{rcp_arg} - {cc_ref_year})")
            tc_clim = tc.apply_climate_scenario_knu(ref_year=int(cc_ref_year), rcp_scenario=rcp_arg)
            out = f"tropcyc/{basin}/TCs_{basin}_{year}_{climate_scenario}_split_{split}.hdf5"
            logger.info(f"Writing to {out}")
            tc_clim.write_hdf5(out)
