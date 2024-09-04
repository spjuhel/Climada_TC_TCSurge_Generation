import sys
import logging, traceback
import os
import copy

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
outdir = Path(snakemake.output[0])
max_tracks = snakemake.params.max_tracks
basin = snakemake.wildcards.genesis_basin
year = snakemake.wildcards.tracks_year

logger.info(f"Loading TC tracks from {snakemake.input.tracks}")

if os.stat(snakemake.input[0]).st_size == 0:
    logger.info(
        f"File is empty, which probably means there is no track data for this basin-year. Ignoring."
    )
    Path(snakemake.output[0]).mkdir(parents=True, exist_ok=True)
else:
    tracks = TCTracks.from_hdf5(snakemake.input.tracks)
    logger.info(f"There are {len(tracks.data)} tracks.")
    Path(snakemake.output[0]).mkdir(parents=True, exist_ok=True)

    split = 1
    for n in range(0, tracks.size, max_tracks):
        logger.info(f"Splitting {n}:{n+max_tracks} tracks")
        tr = copy.deepcopy(tracks)
        tr.data = tr.data[n : n + max_tracks]
        filename = (
            outdir
            / f"IBTracs_{snakemake.config['nsynth']}synth_{basin}_{year}_split_{split}.hdf5"
        )
        logger.info(f"Writing to {filename}")
        tr.write_hdf5(filename)
        split += 1
