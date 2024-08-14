import os
import sys
import logging, traceback
from pathlib import Path

import numpy as np

from climada.hazard import TCTracks, TropCyclone, Centroids
from pathos.pools import ProcessPool as Pool

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

if os.stat(snakemake.input[0]).st_size == 0
    logger.info(f"File is empty, which probably means there is no track data for this basin-year. Ignoring")
    Path(snakemake.output[0]).touch()
else:
    pool = Pool(nodes=snakemake.threads)
    logger.info(f"Created multiprocess pool with {pool.ncpus} cpus.")
    tracks = TCTracks.from_hdf5(snakemake.input.tracks)
    logger.info(f"There are {len(tracks.data)} tracks.")

    logger.info(f"Loading global centroids from {snakemake.input.centroids}")
    cent = Centroids.from_hdf5(snakemake.input.centroids)

    logger.info(f"Selecting centroids extent from tracks with buffer={snakemake.params.buf}")
    cent_tracks = cent.select(extent=tracks.get_extent(snakemake.params.buf))

    logger.info(f"Computing TC wind-fields")
    tc = TropCyclone.from_tracks(tracks, centroids=cent_tracks, pool=pool, max_memory_gb=snakemake.params.max_memory_gb)
    freq_corr = 1/snakemake.config["nsynth"]
    tc.frequency = np.ones(tc.event_id.size)*freq_corr
    pool.close()
    pool.join()

    logger.info(f"Writing to {snakemake.output[0]}")
    tc.write_hdf5(snakemake.output[0])
