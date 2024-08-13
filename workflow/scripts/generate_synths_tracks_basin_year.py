import sys
import logging, traceback
from pathlib import Path

from climada.hazard import TCTracks

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
logger.info(f"Generating Synth TC tracks for genesis basin {snakemake.wildcards.basin} for {snakemake.wildcards.year}")

logger.info(f"Loading TC tracks from {snakemake.input[0]}")
tracks = TCTracks.from_hdf5(snakemake.input[0])

if not tracks.data:
    logger.info(f"No tracks found for this period. Returning empty file")
    Path(snakemake.output[0]).touch()

else:
    tracks.calc_perturbed_trajectories(nb_synth_tracks=snakemake.params.nsynth)

    logger.info(f"Writing to {snakemake.output[0]}")
    tracks.write_hdf5(snakemake.output[0])
