import sys, os
import logging, traceback

from climada.hazard import Hazard

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
logger.info(f"Regrouping")
tcs = Hazard.concat([Hazard.from_hdf5(tcfile) for tcfile in snakemake.input if os.stat(tcfile).st_size!=0])
logger.info(f"Writing to {snakemake.output[0]}")
tcs.write_hdf5(snakemake.output[0])
