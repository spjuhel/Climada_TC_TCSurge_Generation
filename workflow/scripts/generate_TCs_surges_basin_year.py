import sys
import logging, traceback
from pathlib import Path

import pandas as pd
import geopandas as gpd
import xarray as xr

from climada.hazard import TropCyclone
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


basin = snakemake.wildcards.genesis_basin
climate_scenario = snakemake.wildcards.climate_scenario
split = snakemake.wildcards.split
ssp = snakemake.wildcards.ssp
tropcyc = snakemake.input.tropcyc
slr = snakemake.input.slr[0]
slr_year = int(snakemake.wildcards.slr_year)
dem_topo_path = snakemake.input.dem

# Install exception handler
sys.excepthook = handle_exception
logger.info(f"Computing TC Surge events for genesis basin {basin}")

logger.info(f"Reading TC events from {tropcyc}")
tc = TropCyclone.from_hdf5(tropcyc)

if spp=="nossp":
    assert slr_year=="no"
    logger.info(f"Will use DEM data from {dem_topo_path}")
    logger.info(f"Computing surges from TCs")
    ts_rescaled_slr = TCSurgeBathtub.from_tc_winds(
        tc, dem_topo_path, higher_res=snakemake.params.higher_res
    )
    logger.info(f"Writing to {snakemake.output[0]}")
    ts_rescaled_slr.write_hdf5(snakemake.output[0])
else:
    logger.info(f"Reading SLR dataframe from {slr} for year {slr_year}")
    slr_data = xr.open_dataset(slr)
    df = slr_data.to_dataframe()
    slr_data.close()

    df = df.loc[pd.IndexSlice[:, :, slr_year]]
    df.reset_index(drop=True, inplace=True)
    df = df[~df["sea_level_change"].isna()]
    gdf = gpd.GeoDataFrame(
        data=df, geometry=gpd.points_from_xy(df["lon"], df["lat"]), crs="EPSG:4326"
    )
    gdf = gdf.fillna(0)

    logger.info(f"Will use DEM data from {dem_topo_path}")
    ## read by the method nothing to do here

    logger.info(f"Computing surges from TCs")
    ts_rescaled_slr = TCSurgeBathtub.from_tc_winds(
        tc, dem_topo_path, higher_res=snakemake.params.higher_res, sea_level_rise_gdf=gdf
    )

    logger.info(f"Writing to {snakemake.output[0]}")
    ts_rescaled_slr.write_hdf5(snakemake.output[0])
