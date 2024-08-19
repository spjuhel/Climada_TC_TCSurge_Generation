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

logger.info(f"Getting TC tracks for genesis basin {snakemake.wildcards.basin} for {snakemake.wildcards.year}")

tracks = TCTracks.from_ibtracs_netcdf(year_range=(int(snakemake.wildcards.year), int(snakemake.wildcards.year)), genesis_basin=snakemake.wildcards.basin, estimate_missing=True)

if not tracks.data:
    logger.info(f"No tracks found for this period. Returning empty file")
    Path(snakemake.output[0]).touch()

else:
    logger.info(f"Interpolating tracks to {snakemake.params.timestep} hours steps")
    tracks.equal_timestep(snakemake.params.timestep)

    logger.info(f"Writing to {snakemake.output[0]}")
    tracks.write_hdf5(snakemake.output[0])
