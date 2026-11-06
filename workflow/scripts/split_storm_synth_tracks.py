import copy
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

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.StreamHandler())

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    LOGGER.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )

# Install exception handler
sys.excepthook = handle_exception

def main(
    tracks_path,
    output_dir_path,
    output_name_format,
    max_tracks
):
    """Compute and write the TC hazard for one STORM basin/member.

    Parameters
    ----------
    tracks_path : str or pathlib.Path
        STORM track file.
    output_dir_path : str or pathlib.Path
        Where to write the resulting splitted Tracks (hdf5).
    """
    LOGGER.info("Loading STORM tracks from %s", tracks_path)
    tracks = TCTracks.from_hdf5(tracks_path)
    if not tracks.data:
        LOGGER.info("No tracks found for this file. Ignoring.")
        Path(output_dir_path).touch()

    LOGGER.info(f"There are {len(tracks.data)} tracks.")
    Path(output_dir_path).mkdir(parents=True, exist_ok=True)

    split = 1
    for n in range(0, tracks.size, max_tracks):
        LOGGER.info(f"Splitting {n}:{n+max_tracks} tracks")
        tr = copy.deepcopy(tracks)
        tr.data = tr.data[n : n + max_tracks]
        filename = (
            Path(output_dir_path)
            / output_name_format.format(split)
        )
        LOGGER.info(f"Writing to {filename}")
        tr.write_hdf5(filename)
        split += 1

if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    main(
        tracks_path=snakemake.input.tracks,
        output_dir_path=snakemake.output[0],
        output_name_format=snakemake.params.output_name_format,
        max_tracks=snakemake.params.max_tracks,
    )
