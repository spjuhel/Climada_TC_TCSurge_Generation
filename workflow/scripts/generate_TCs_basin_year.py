import sys
import logging, traceback
from pathlib import Path

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
tracks = TCTracks.from_hdf5(snakemake.input.tracks)

if not tracks.data:
    logger.info(f"No tracks found for this period. Returning empty file")
    Path(snakemake.output[0]).touch()

else:
    logger.info(f"Loading global centroids from {snakemake.input.global_cent}")
    cent = Centroids.from_hdf5(snakemake.input.global_cent)

    logger.info(f"Selecting centroids extent from tracks with buffer={snakemake.params.buf}")
    cent_tracks = cent.select(extent=tracks.get_extent(snakemake.params.buf))

    logger.info(f"Computing TC wind-fields")
    tc = TropCyclone.from_tracks(tracks, centroids=cent_tracks)

    logger.info(f"Writing to {snakemake.output[0]}")
    tc.write_hdf5(snakemake.output[0])