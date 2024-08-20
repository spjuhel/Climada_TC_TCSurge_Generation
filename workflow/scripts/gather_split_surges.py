import sys
import logging, traceback

from climada_petals.hazard import TCSurgeBathtub

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
tcs = TCSurgeBathtub.concat([TCSurgeBathtub.from_hdf5(tcfile) for tcfile in snakemake.input])
logger.info(f"Writing to {snakemake.output[0]}")
tcs.write_hdf5(snakemake.output[0])
