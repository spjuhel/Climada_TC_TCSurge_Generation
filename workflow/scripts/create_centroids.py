import sys
import logging, traceback

from climada.hazard import Centroids

def create_litpop_matching_centroids(res=0.041666659999975764):
    centroids_land = Centroids.from_pnt_bounds(
        (-180 + res / 2, -90 + res / 2, 180 - res / 2, 90), res
    )
    centroids_land.set_on_land()

    centroids_water = Centroids.from_pnt_bounds((-180, -90, 180, 90), 0.5)
    centroids_water.set_on_land()
    cent = centroids_land.select(sel_cen=centroids_land.on_land).union(
        centroids_water.select(sel_cen=~centroids_water.on_land)
    )
    return cent

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


logger.info(f"Creating global centroids ")
cent = create_litpop_matching_centroids()
logger.info(f"Writing to {snakemake.output[0]}")
cent.write_hdf5(snakemake.output[0])
