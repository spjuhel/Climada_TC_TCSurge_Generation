import sys
import logging, traceback

from climada.hazard import TCTracks, Centroids, TropCyclone

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

logger.info(f"Getting TC tracks for genesis basin {snakemake.wildcards.basin} for {snakemake.wildcards.years}")

tracks = TCTracks.from_ibtracs_netcdf(year_range=(int(snakemake.wildcards.start), int(snakemake.wildcards.end)), genesis_basin=snakemake.wildcards.basin)

logger.info(f"Loading global centroids from {snakemake.input.global_cent}")
cent = Centroids.from_hdf5(snakemake.input.global_cent)

logger.info(f"Selecting centroids extent from tracks with buffer={snakemake.params.buf}")
cent_tracks = cent.select(extent=tracks.get_extent(snakemake.params.buf))

logger.info(f"Interpolating tracks to {snakemake.params.timestep} hours steps")
tracks.equal_timestep(snakemake.params.timestep)

logger.info(f"Computing TC wind-fields")
tc = TropCyclone.from_tracks(tracks, centroids=cent_tracks)

logger.info(f"Writing to {snakemake.output[0]}")
tc.write_hdf5(snakemake.output[0])
