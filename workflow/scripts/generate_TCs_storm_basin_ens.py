"""Compute tropical cyclone wind fields from STORM tracks.

Snakemake script. For one basin and one ensemble member of the STORM
present-climate simulations, load the tracks, interpolate them to a fixed
time step and compute the wind fields on the global centroids.

Adapted from a script by Simona Meiler (2024-02-26).
"""

import logging
from pathlib import Path

from climada.hazard import Centroids, TCTracks, TropCyclone

LOGGER = logging.getLogger(__name__)


def load_storm_tracks(tracks_path, time_step_h=1.0):
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
    tracks.equal_timestep(time_step_h=time_step_h)
    return tracks


def main(
    tracks_path,
    centroids_path,
    output_path,
    wind_model="H08",
    time_step_h=1.0,
    max_memory_gb=14,
):
    """Compute and write the TC hazard for one STORM basin/member.

    Parameters
    ----------
    tracks_path : str or pathlib.Path
        STORM track file.
    centroids_path : str or pathlib.Path
        CLIMADA centroids file (hdf5).
    output_path : str or pathlib.Path
        Where to write the resulting `TropCyclone` hazard (hdf5).
    wind_model : str, default "H08"
        Wind field model passed to `TropCyclone.from_tracks`
        (e.g. "H08" for Holland 2008, "ER11" for Emanuel and Rotunno 2011).
    time_step_h : float, default 1.0
        Time step in hours of the interpolated tracks.
    max_memory_gb : float, default 14
        Memory budget passed to `TropCyclone.from_tracks`.
    """
    LOGGER.info("Loading STORM tracks from %s", tracks_path)
    tracks = load_storm_tracks(tracks_path, time_step_h=time_step_h)

    LOGGER.info("Loading centroids from %s", centroids_path)
    centroids = Centroids.from_hdf5(str(centroids_path))

    LOGGER.info("Computing wind fields with model %s", wind_model)
    tc_hazard = TropCyclone.from_tracks(
        tracks,
        centroids=centroids,
        model=wind_model,
        max_memory_gb=max_memory_gb,
    )
    tc_hazard.check()

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
        tracks_path=snakemake.input.tracks,
        centroids_path=snakemake.input.centroids,
        output_path=snakemake.output[0],
        wind_model=snakemake.params.wind_model,
        time_step_h=snakemake.params.time_step_h,
        max_memory_gb=snakemake.params.max_memory_gb,
    )
