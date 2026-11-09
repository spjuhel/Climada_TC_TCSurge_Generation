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


def main(
    tracks_path,
    centroids_path,
    output_path,
    buffer,
    wind_model="H08",
    temp_res=.1,
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
    tracks = TCTracks.from_hdf5(tracks_path)
    LOGGER.info(f"There are {len(tracks.data)} tracks.")

    LOGGER.info("Loading centroids from %s", centroids_path)
    centroids = Centroids.from_hdf5(str(centroids_path))

    LOGGER.info(f"Selecting centroids extent from tracks with buffer={buffer}")
    cent_tracks = centroids.select(extent=tracks.get_extent(buffer))

    LOGGER.info(f"Increasing temporal resolution of tracks to {temp_res}")
    tracks.equal_timestep(temp_res)

    LOGGER.info("Computing wind fields with model %s", wind_model)
    tc_hazard = TropCyclone.from_tracks(
        tracks,
        centroids=cent_tracks,
        model=wind_model,
        max_memory_gb=max_memory_gb,
    )

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
        buffer=snakemake.params.buffer,
        wind_model=snakemake.params.wind_model,
        temp_res=snakemake.params.time_step_h,
        max_memory_gb=(snakemake.resources.get("mem_mb_per_cpu")-2000)/10,
    )
