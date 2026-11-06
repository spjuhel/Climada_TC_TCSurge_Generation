"""Compute tropical cyclone wind fields from STORM tracks.

Snakemake script. For one basin of the STORM present-climate simulations,
contenate all ensemble members windfield and sets the frequency.

Adapted from scripts by Simona Meiler (2024-02-26).
"""

import logging
from pathlib import Path

from climada.hazard import TropCyclone

LOGGER = logging.getLogger(__name__)

def concatenate_members(files):
    trop_cyc = TropCyclone.from_hdf5(files[0])
    for file in files[1:]:
        trop_cyc.append(TropCyclone.from_hdf5(file))

    return trop_cyc

def main(files, output_path):
    """."""
    LOGGER.info("Loading members from %s", files.parent)
    tc_hazard = concatenate_members(files)
    LOGGER.info("Writing hazard to %s", output_path)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    tc_hazard.write_hdf5(str(output_path))


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    main(
        files=snakemake.input.members,
        output_path=snakemake.output[0],
    )
