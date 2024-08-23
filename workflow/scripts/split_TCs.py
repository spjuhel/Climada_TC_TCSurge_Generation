import sys
import logging, traceback
import os
import copy

from pathlib import Path

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

tropcyc = snakemake.input.tropcyc
logger.info(f"Loading TCs from {tropcyc}")

outdir = Path(snakemake.output[0])
max_tcs = snakemake.params.max_tracks

if os.stat(snakemake.input[0]).st_size == 0:
    logger.info(
        f"File is empty, which probably means there is no track data for this basin-year. Ignoring."
    )
    Path(snakemake.output[0]).mkdir(parents=True, exist_ok=True)
else:
    tcs = TropCyclone.from_hdf5(tropcyc)
    logger.info(f"There are {tropcyc.size} TCs.")
    Path(snakemake.output[0]).mkdir(parents=True, exist_ok=True)

    split = 1
    for n in range(0, tcs.size, max_tcs):
        logger.info(f"Splitting {n}:{n+max_tcs} TCs")
        tc = tcs.select(event_id=[n+1 : n + max_tcs+1])
        filename = (
            outdir
            / f"TCs_{snakemake.wildcards.basin}_{snakemake.wildcards.climate_scenario}_split_{split}.hdf5"
        )
        logger.info(f"Writing to {filename}")
        tc.write_hdf5(filename)
        split += 1
