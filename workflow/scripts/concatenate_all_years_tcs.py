import os
import sys
import logging, traceback

from climada.hazard import TropCyclone

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

logger.info(f"Regrouping all periods for basin: {basin}")
tcfiles = [fname for fname in snakemake.input if os.stat(fname).st_size != 0]
tcs = TropCyclone.concat([TropCyclone.from_hdf5(tcfile) for tcfile in tcfiles])
tcs.frequency /= (snakemake.config["end"] - snakemake.config["start"])
logger.info(f"Writing to {snakemake.output[0]}")
tcs.write_hdf5(snakemake.output[0])
