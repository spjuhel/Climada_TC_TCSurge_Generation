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


def load_storm_tracks(tracks_path, timestep=1.0):
    """Load STORM tracks and interpolate them to a fixed time step.

    Parameters
    ----------
    tracks_path : str or pathlib.Path
        Path to a STORM text file (``STORM_DATA_IBTRACS_*.txt``).
    time_step_h : float, default 1.0
        Time step in hours of the interpolated tracks.

    Returns
    -------
    climada.hazard.TCTracks
        Tracks with an equal time step.
    """
    tracks = TCTracks.from_simulations_storm(str(tracks_path))
    tracks.equal_timestep(time_step_h=timestep)
    if tracks.data:
        LOGGER.info(f"Interpolating tracks to {timestep} hours steps")
        tracks.equal_timestep(timestep)

    return tracks

def main(
    tracks_path,
    output_path,
    time_step_h=1.0,
):
    """Compute and write the TC hazard for one STORM basin/member.

    Parameters
    ----------
    tracks_path : str or pathlib.Path
        STORM track file.
    output_path : str or pathlib.Path
        Where to write the resulting Tracks (hdf5).
    time_step_h : float, default 1.0
        Time step in hours of the interpolated tracks.
    """
    LOGGER.info("Loading STORM tracks from %s", tracks_path)
    tracks = load_storm_tracks(tracks_path, timestep=time_step_h)
    if not tracks.data:
        LOGGER.info("No tracks found for this period. Returning empty file")
        Path(output_path).touch()
    LOGGER.info("Writing STORM tracks to %s", output_path)
    tracks.write_hdf5(output_path)


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    main(
        tracks_path=snakemake.input.tracks,
        output_path=snakemake.output[0],
        time_step_h=snakemake.params.timestep,
    )
