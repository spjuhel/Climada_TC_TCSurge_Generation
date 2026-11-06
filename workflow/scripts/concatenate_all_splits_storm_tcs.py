import logging
import os

from climada.hazard import TropCyclone

LOGGER = logging.getLogger(__name__)


def main(
    tcsfiles,
    output_path,
):
    LOGGER.info("Loading STORM tropical cyclones from %s", tcsfiles)
    tcfiles = [fname for fname in tcsfiles if os.stat(fname).st_size != 0]
    tcs_set = [TropCyclone.from_hdf5(tcfile) for tcfile in tcfiles]
    tcs_set = [tc for tc in tcs_set if tc.size>0]
    tcs = TropCyclone.concat(tcs_set)
    LOGGER.info(f"Writing to {output_path}")
    tcs.write_hdf5(output_path)


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    main(
        tcsfiles=snakemake.input.tracks,
        output_path=snakemake.output[0],
    )
