import os
import sys
import logging, traceback
from pathlib import Path

import numpy as np
import copy

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
logger.info(f"Computing TC events for genesis basin {snakemake.wildcards.basin} for {snakemake.wildcards.year}")

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
    tclist = []
    for n in range(0, tracks.size, snakemake.params.batch_size):
        tr = copy.deepcopy(tracks)
        tr.data = tr.data[n:n+snakemake.params.batch_size]
        tr.equal_timestep(0.5)
        tc = TropCyclone.from_tracks(tr, centroids=cent_tracks, max_memory_gb=snakemake.params.max_memory_gb)
        freq_corr = 1/snakemake.config["nsynth"]
        tc.frequency = np.ones(tc.event_id.size)*freq_corr
        tclist.append(tc)

    tc = TropCyclone.concat(tclist)
    logger.info(f"Writing to {snakemake.output[0]}")
    tc.write_hdf5(snakemake.output[0])
